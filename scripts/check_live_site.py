#!/usr/bin/env python3
"""Check the live site: every internal page, figure, style sheet, download and #anchor (external links too).

Crawls every HTML page reachable from URL under the same scheme, host and path prefix, fetches every href
and src that resolves there (GET, 200 expected) and checks that each #fragment used by a link exists as an
id (or <a name>) in its target page. --external also fetches every other http(s) link with a desktop browser
User-Agent (redirects followed, 20-s timeout, one retry): 2xx/3xx is ok; 401, 403 or 429 from a known
publisher or DOI resolver (PUBLISHERS) is "blocked to scripts" (check it in a browser); anything else fails.
Standard library only. Exit status 1 on an internal failure or a failed external link, not on a blocked one.
After a push, GitHub's CDN may serve the previous build for a few minutes.

Usage: .venv/bin/python scripts/check_live_site.py https://khan-laiba.github.io/opm-squid-meg/ [--external]
"""
from __future__ import annotations

import argparse
import http.client
import http.cookiejar
import re
import sys
import time
import urllib.error
import urllib.request
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
from html.parser import HTMLParser
from typing import NamedTuple
from urllib.parse import quote, unquote, urldefrag, urljoin, urlsplit, urlunsplit

TIMEOUT = 20  # s per request
RETRY_DELAY = 2.0  # s before the one retry of a failed request
WORKERS = 8  # parallel requests
SCRIPT_HEADERS = {"User-Agent": "check_live_site.py (Python urllib)"}
BROWSER_HEADERS = {"User-Agent": "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 (KHTML, like Gecko) "
                                 "Chrome/129.0.0.0 Safari/537.36",
                   "Accept": "text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8",
                   "Accept-Language": "en-US,en;q=0.9"}
# DOI and handle resolvers, publishers, preprint servers and literature databases (with their subdomains) that
# often answer scripts with 401, 403 or 429 although the link works in a browser
PUBLISHERS = ("doi.org", "handle.net", "biorxiv.org", "medrxiv.org", "arxiv.org", "researchsquare.com", "ncbi.nlm.nih.gov",
              "europepmc.org", "wiley.com", "sciencedirect.com", "elsevier.com", "cell.com", "thelancet.com", "springer.com",
              "nature.com", "biomedcentral.com", "iop.org", "tandfonline.com", "sagepub.com", "oup.com", "science.org",
              "pnas.org", "frontiersin.org", "mdpi.com", "plos.org", "ieee.org", "cambridge.org", "lww.com", "neurology.org",
              "jneurosci.org", "physiology.org", "jamanetwork.com", "nejm.org", "bmj.com", "karger.com", "thieme-connect.com",
              "direct.mit.edu", "elifesciences.org", "aip.org", "aps.org", "acs.org", "optica.org", "annualreviews.org",
              "royalsocietypublishing.org", "degruyter.com", "hindawi.com", "jstor.org")


class Result(NamedTuple):
    status: int | None  # final HTTP status; None if there was no response
    url: str  # final URL, after redirects
    ctype: str  # content type of a successful response
    text: str | None  # the page, when an HTML response was read as one
    error: str  # reason of a failure


class _Page(HTMLParser):
    """The href and src values, the ids (and <a name>) and the <base href> of one HTML page."""

    def __init__(self):
        super().__init__()
        self.refs, self.ids, self.base = [], set(), None

    def handle_starttag(self, tag, attrs):
        a = {k: (v or "").strip() for k, v in attrs}
        if a.get("id"):
            self.ids.add(a["id"])
        if tag == "a" and a.get("name"):
            self.ids.add(a["name"])
        if tag == "base":
            self.base = self.base or a.get("href") or None
            return
        if tag == "link" and {"preconnect", "dns-prefetch"} & set(a.get("rel", "").lower().split()):
            return  # a hint naming an origin, not a resource
        self.refs += [a[k] for k in ("href", "src") if a.get(k)]


def norm(url: str) -> str:
    """``url`` without its fragment, scheme and host in lower case, unsafe characters percent-encoded."""
    s = urlsplit(urldefrag(url)[0])
    return urlunsplit((s.scheme.lower(), s.netloc.lower(), quote(s.path or "/", safe="/%!$&'()*+,;=:@~"),
                       quote(s.query, safe="/?%!$&'()*+,;=:@~"), ""))


def fetch(url: str, headers: dict, read: bool = False, cookies: bool = False) -> Result:
    """GET ``url``, following redirects; with ``read`` the whole body is downloaded (a truncated one fails) and an
    HTML page is kept as text."""
    opener = urllib.request.build_opener(*([urllib.request.HTTPCookieProcessor(http.cookiejar.CookieJar())] if cookies else []))
    try:
        with opener.open(urllib.request.Request(url, headers=headers), timeout=TIMEOUT) as r:
            ctype, body = r.headers.get_content_type(), r.read() if read else None
            text = body.decode(r.headers.get_content_charset() or "utf-8", "replace") if read and ctype == "text/html" else None
            return Result(r.status, r.geturl(), ctype, text, "")
    except urllib.error.HTTPError as e:
        e.close()
        return Result(e.code, e.filename or url, "", None, str(e.reason))
    except (urllib.error.URLError, http.client.HTTPException, OSError, ValueError) as e:
        return Result(None, url, "", None, str(getattr(e, "reason", e)))


def get(url: str, headers: dict, **kw) -> Result:
    """fetch() with one retry after a failure (no response, or a 4xx or 5xx status)."""
    r = fetch(url, headers, **kw)
    if r.status is None or r.status >= 400:
        time.sleep(RETRY_DELAY)
        r = fetch(url, headers, **kw)
    return r


def is_publisher(url: str) -> bool:
    host = (urlsplit(url).hostname or "").lower()
    return any(host == d or host.endswith("." + d) for d in PUBLISHERS)


def classify(status: int | None, url: str, final: str) -> str:
    """'ok' (2xx or 3xx), 'blocked' (401, 403 or 429 from a publisher or DOI resolver, before or after redirects)
    or 'failed'."""
    if status is not None and 200 <= status < 400:
        return "ok"
    if status in (401, 403, 429) and (is_publisher(url) or is_publisher(final)):
        return "blocked"
    return "failed"


def check(url: str, external: bool = False) -> dict:
    """Crawl the site at ``url`` and check its links (the external ones only if ``external``).

    Returns root (the site's prefix), pages (distinct HTML pages; dir/ and dir/index.html count once), assets
    (the other internal URLs), anchors (target and fragment pairs checked), failures [(message, linking pages)]
    and external {url: verdict (None if not checked), status, final, error, linked_from}."""
    url = url if "://" in url else "https://" + url
    first = get(norm(url), SCRIPT_HEADERS, read=True)
    start = norm(first.url if first.status == 200 else url)  # after redirects, e.g. .../repository -> .../repository/
    s = urlsplit(start)
    root = urlunsplit((s.scheme, s.netloc, s.path[: s.path.rfind("/") + 1], "", ""))
    results, pages, failures = {start: first}, {}, []
    linked, anchors, ext = defaultdict(set), defaultdict(set), defaultdict(set)  # target -> pages linking to it

    def visit(page_url, r):
        page = _Page()
        page.feed(r.text)
        pages[page_url] = page
        base = urljoin(r.url, page.base) if page.base else r.url
        for ref in page.refs:
            try:
                target, frag = urldefrag(urljoin(base, ref))
                if urlsplit(target).scheme not in ("http", "https"):
                    continue  # mailto:, javascript:, data: ...
                target = norm(target)
            except ValueError:
                failures.append((f"malformed link {ref!r}", {page_url}))
                continue
            if target.startswith(root):
                linked[target].add(page_url)
                if frag:
                    anchors[target, frag].add(page_url)
            else:
                ext[target].add(page_url)

    if first.status == 200 and first.ctype == "text/html":
        visit(start, first)
    elif first.status == 200:
        failures.append((f"not an HTML page: {start}", set()))
    while todo := sorted(set(linked) - set(results)):
        with ThreadPoolExecutor(WORKERS) as pool:
            got = list(pool.map(lambda u: get(u, SCRIPT_HEADERS, read=True), todo))
        for u, r in zip(todo, got, strict=True):
            results[u] = r
            if r.status == 200 and r.ctype == "text/html" and norm(r.url).startswith(root):
                visit(u, r)
    for u, r in sorted(results.items()):
        if r.status != 200:
            failures.append((f"{r.status or 'no response'} {u}" + (f" ({r.error})" if r.error else ""), linked.get(u, set())))
    n_anchors = 0
    for (u, frag), refs in sorted(anchors.items()):
        if u not in pages:
            continue  # a failed target (reported above) or not an HTML page
        n_anchors += 1
        if frag.lower() != "top" and not {frag, unquote(frag)} & pages[u].ids:
            failures.append((f"missing anchor {u}#{frag}", refs))
    checked = {}
    if external:
        urls = sorted(ext)
        with ThreadPoolExecutor(WORKERS) as pool:
            checked = dict(zip(urls, pool.map(lambda u: get(u, BROWSER_HEADERS, cookies=True), urls), strict=True))
    out = {}
    for u, refs in sorted(ext.items()):
        r = checked.get(u)
        out[u] = dict(verdict=classify(r.status, u, r.url) if r else None, status=r and r.status, final=r and r.url,
                      error=r and r.error, linked_from=sorted(refs))
    return dict(root=root, pages=len({re.sub(r"/index\.html?$", "/", u) for u in pages}), assets=len(results) - len(pages),
                anchors=n_anchors, failures=failures, external=out)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("url", help="address of the site, e.g. https://<owner>.github.io/<repository>/")
    ap.add_argument("--external", action="store_true", help="also check every external http(s) link")
    args = ap.parse_args(argv)
    rep = check(args.url, external=args.external)
    root, ext = rep["root"], rep["external"]

    def where(refs):
        return ", ".join(sorted((r[len(root):] or "./") if r.startswith(root) else r for r in refs)) or "the start URL"

    for msg, refs in rep["failures"]:
        print(f"FAILED   {msg}   <- {where(refs)}")
    for u, v in ext.items():
        if v["verdict"] in ("blocked", "failed"):
            tag = "BLOCKED to scripts (check in a browser)" if v["verdict"] == "blocked" else "FAILED  "
            why = f"{v['status']} {v['error']}" if v["status"] else f"no response ({v['error']})"
            to = f" -> {v['final']}" if v["final"] and norm(v["final"]) != u else ""
            print(f"{tag} {why}: {u}{to}   <- {where(v['linked_from'])}")
    n = Counter(v["verdict"] for v in ext.values())
    print(f"{root}: {rep['pages']} pages, {rep['assets']} internal assets, {rep['anchors']} anchors checked, "
          f"{len(rep['failures'])} internal failures")
    print(f"external links: {len(ext)}" + (f"; {n['ok']} ok, {n['blocked']} blocked to scripts (check in a browser), "
                                          f"{n['failed']} failed" if args.external else ", not checked (use --external)"))
    return 1 if rep["failures"] or n["failed"] else 0


if __name__ == "__main__":
    sys.exit(main())
