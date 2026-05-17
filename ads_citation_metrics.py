#!/usr/bin/env python3
"""Compute ADS citation metrics for all papers in a library.

Outputs:
1) Total citations across all papers in the library.
2) Total citations excluding self-citations.

Self-citation rule used here:
- A citation is treated as self-citation when the citing paper and cited paper
  share at least one normalized author key.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
import sys
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import dataclass
from typing import Any, Dict, Iterable, List, Optional, Set


ADS_API_BASE = "https://api.adsabs.harvard.edu/v1"
DEFAULT_ROWS = 2000
CHUNK_SIZE = 100


def chunks(items: List[str], size: int) -> Iterable[List[str]]:
    for i in range(0, len(items), size):
        yield items[i : i + size]


def extract_library_id(value: str) -> str:
    value = value.strip()
    if not value:
        raise ValueError("Library URL or ID is empty.")
    if "/public-libraries/" in value:
        return value.rstrip("/").split("/")[-1]
    return value


def normalize_author(author: str) -> str:
    """Normalize an ADS author string to a comparable key."""
    author = author.strip().lower()
    if not author:
        return ""

    # Keep alphanumerics, comma, and spaces for robust parsing.
    cleaned = re.sub(r"[^a-z0-9, ]+", " ", author)
    cleaned = re.sub(r"\s+", " ", cleaned).strip()

    if "," in cleaned:
        last, first = cleaned.split(",", 1)
        last = last.strip()
        first_tokens = first.strip().split()
        first_initial = first_tokens[0][0] if first_tokens and first_tokens[0] else ""
        return f"{last}|{first_initial}"

    tokens = cleaned.split()
    if len(tokens) == 1:
        return f"{tokens[0]}|"
    last = tokens[-1]
    first_initial = tokens[0][0] if tokens[0] else ""
    return f"{last}|{first_initial}"


def normalized_author_set(authors: List[str]) -> Set[str]:
    return {key for key in (normalize_author(a) for a in authors) if key}


def parse_publication_year(doc: Dict) -> Optional[int]:
    year = doc.get("year")
    if isinstance(year, int):
        return year
    if isinstance(year, str) and year.isdigit():
        return int(year)

    pubdate = doc.get("pubdate")
    if isinstance(pubdate, str):
        match = re.search(r"\b(\d{4})\b", pubdate)
        if match:
            return int(match.group(1))
    return None


def parse_title(doc: Dict) -> str:
    title = doc.get("title")
    if isinstance(title, list) and title:
        first = title[0]
        if isinstance(first, str):
            return first
    if isinstance(title, str):
        return title
    return ""


@dataclass
class PaperRecord:
    bibcode: str
    title: str
    authors: Set[str]
    citation_count: int
    citation_bibcodes: List[str]
    publication_year: Optional[int]


class AdsClient:
    def __init__(self, token: str) -> None:
        self.token = token

    def _get_json(self, path: str, query: Dict[str, str] | None = None) -> Dict:
        url = f"{ADS_API_BASE}{path}"
        if query:
            url = f"{url}?{urllib.parse.urlencode(query)}"

        req = urllib.request.Request(url)
        req.add_header("Authorization", f"Bearer {self.token}")
        req.add_header("Accept", "application/json")
        req.add_header(
            "User-Agent",
            "ads-citation-metrics/1.0 (+https://ui.adsabs.harvard.edu/)",
        )

        body = ""
        try:
            with urllib.request.urlopen(req, timeout=60) as resp:
                body = resp.read().decode("utf-8")
        except urllib.error.HTTPError as exc:
            details = exc.read().decode("utf-8", errors="replace")
            # Some environments get blocked by ADS WAF when using urllib,
            # but succeed with curl and the same token.
            if exc.code == 405 or "Human Verification" in details:
                return self._get_json_via_curl(url)
            raise RuntimeError(f"ADS HTTP {exc.code}: {details}") from exc
        except urllib.error.URLError as exc:
            raise RuntimeError(f"Network error reaching ADS API: {exc}") from exc

        try:
            return json.loads(body)
        except json.JSONDecodeError as exc:
            if "Human Verification" in body:
                return self._get_json_via_curl(url)
            raise RuntimeError("ADS API returned non-JSON response.") from exc

    def _get_json_via_curl(self, url: str) -> Dict:
        command = [
            "curl",
            "-sS",
            "--fail-with-body",
            url,
            "-H",
            f"Authorization: Bearer {self.token}",
            "-H",
            "Accept: application/json",
            "-H",
            "User-Agent: ads-citation-metrics/1.0 (+https://ui.adsabs.harvard.edu/)",
        ]
        try:
            completed = subprocess.run(
                command,
                check=True,
                capture_output=True,
                text=True,
            )
        except FileNotFoundError as exc:
            raise RuntimeError(
                "curl is required for fallback API access but was not found on PATH."
            ) from exc
        except subprocess.CalledProcessError as exc:
            stderr = (exc.stderr or "").strip()
            stdout = (exc.stdout or "").strip()
            details = stderr or stdout or "Unknown curl error."
            raise RuntimeError(f"ADS request failed via curl: {details}") from exc

        body = completed.stdout
        try:
            return json.loads(body)
        except json.JSONDecodeError as exc:
            raise RuntimeError("ADS API returned non-JSON response via curl.") from exc

    def get_library_bibcodes(self, library_id: str, rows: int = DEFAULT_ROWS) -> List[str]:
        data = self._get_json(
            f"/biblib/libraries/{library_id}",
            query={"rows": str(rows)},
        )
        docs = data.get("documents", [])
        if not isinstance(docs, list):
            raise RuntimeError("Unexpected ADS library response format (documents missing).")
        return [d for d in docs if isinstance(d, str)]

    def fetch_paper_docs(self, bibcodes: List[str], fields: List[str]) -> List[Dict]:
        all_docs: List[Dict] = []
        fl = ",".join(fields)

        for batch in chunks(bibcodes, CHUNK_SIZE):
            escaped = []
            for code in batch:
                code = code.replace("\\", "\\\\").replace('"', '\\"')
                escaped.append(f'"{code}"')
            q = f"bibcode:({ ' OR '.join(escaped) })"

            payload = self._get_json(
                "/search/query",
                query={"q": q, "fl": fl, "rows": str(len(batch))},
            )
            docs = payload.get("response", {}).get("docs", [])
            if isinstance(docs, list):
                all_docs.extend(docs)
        return all_docs


def build_paper_records(client: AdsClient, bibcodes: List[str]) -> Dict[str, PaperRecord]:
    fields = ["bibcode", "title", "author", "citation_count", "citation", "year", "pubdate"]
    docs = client.fetch_paper_docs(bibcodes, fields)

    records: Dict[str, PaperRecord] = {}
    for doc in docs:
        bibcode = doc.get("bibcode")
        if not isinstance(bibcode, str):
            continue

        authors = doc.get("author") or []
        citation_bibcodes = doc.get("citation") or []
        citation_count = doc.get("citation_count") or 0

        if not isinstance(authors, list):
            authors = []
        if not isinstance(citation_bibcodes, list):
            citation_bibcodes = []
        if not isinstance(citation_count, int):
            try:
                citation_count = int(citation_count)
            except (TypeError, ValueError):
                citation_count = 0

        records[bibcode] = PaperRecord(
            bibcode=bibcode,
            title=parse_title(doc),
            authors=normalized_author_set([a for a in authors if isinstance(a, str)]),
            citation_count=citation_count,
            citation_bibcodes=[c for c in citation_bibcodes if isinstance(c, str)],
            publication_year=parse_publication_year(doc),
        )
    return records


def compute_metrics(
    client: AdsClient, library_id: str, year_from: int, year_to: int
) -> Dict[str, Any]:
    library_bibcodes = client.get_library_bibcodes(library_id)
    if not library_bibcodes:
        return {
            "papers_in_library": 0,
            "total_citations": 0,
            "total_citations_excluding_self": 0,
            "top_paper_by_citations_excluding_self": None,
        }

    all_paper_records = build_paper_records(client, library_bibcodes)
    paper_records = {
        bibcode: record
        for bibcode, record in all_paper_records.items()
        if record.publication_year is not None
        and year_from <= record.publication_year <= year_to
    }
    if not paper_records:
        return {
            "papers_in_library": len(library_bibcodes),
            "papers_resolved": len(all_paper_records),
            "papers_in_year_range": 0,
            "year_from": year_from,
            "year_to": year_to,
            "total_citations": 0,
            "self_citations_detected": 0,
            "total_citations_excluding_self": 0,
            "top_paper_by_citations_excluding_self": None,
        }

    total_citations = sum(r.citation_count for r in paper_records.values())

    # Gather all citing bibcodes for author lookup.
    all_citing_bibcodes: Set[str] = set()
    for record in paper_records.values():
        all_citing_bibcodes.update(record.citation_bibcodes)

    citing_author_map: Dict[str, Set[str]] = {}
    if all_citing_bibcodes:
        citing_docs = client.fetch_paper_docs(
            sorted(all_citing_bibcodes),
            fields=["bibcode", "author"],
        )
        for doc in citing_docs:
            bibcode = doc.get("bibcode")
            authors = doc.get("author") or []
            if isinstance(bibcode, str) and isinstance(authors, list):
                citing_author_map[bibcode] = normalized_author_set(
                    [a for a in authors if isinstance(a, str)]
                )

    self_citation_count = 0
    per_paper_citation_stats: Dict[str, Dict[str, Any]] = {}
    for record in paper_records.values():
        self_citing_bibcodes: List[str] = []
        non_self_citing_bibcodes: List[str] = []
        for citing_bibcode in record.citation_bibcodes:
            citing_authors = citing_author_map.get(citing_bibcode)
            if citing_authors and (record.authors & citing_authors):
                self_citation_count += 1
                self_citing_bibcodes.append(citing_bibcode)
            else:
                non_self_citing_bibcodes.append(citing_bibcode)

        per_paper_citation_stats[record.bibcode] = {
            "total_citations": record.citation_count,
            "self_citations": len(self_citing_bibcodes),
            "citations_excluding_self": len(non_self_citing_bibcodes),
            "self_citing_bibcodes": self_citing_bibcodes,
            "non_self_citing_bibcodes": non_self_citing_bibcodes,
        }

    top_record = max(
        paper_records.values(),
        key=lambda r: per_paper_citation_stats[r.bibcode]["citations_excluding_self"],
    )
    top_stats = per_paper_citation_stats[top_record.bibcode]
    top_paper_info = {
        "title": top_record.title or top_record.bibcode,
        "bibcode": top_record.bibcode,
        "publication_year": top_record.publication_year,
        # Explicit counts for quick reading.
        "total_citations": top_stats["total_citations"],
        "self_citations": top_stats["self_citations"],
        "citations_excluding_self": top_stats["citations_excluding_self"],
        "citation_information": top_stats,
    }

    return {
        "papers_in_library": len(library_bibcodes),
        "papers_resolved": len(all_paper_records),
        "papers_in_year_range": len(paper_records),
        "year_from": year_from,
        "year_to": year_to,
        "total_citations": total_citations,
        "self_citations_detected": self_citation_count,
        "total_citations_excluding_self": max(total_citations - self_citation_count, 0),
        "top_paper_by_citations_excluding_self": top_paper_info,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compute ADS citation totals and citation totals excluding self-citations.",
    )
    parser.add_argument(
        "--library",
        required=True,
        help="ADS library ID or full public library URL.",
    )
    parser.add_argument(
        "--token",
        default=os.getenv("ADS_API_TOKEN", ""),
        help="ADS API token. Defaults to ADS_API_TOKEN environment variable.",
    )
    parser.add_argument(
        "--year-from",
        type=int,
        default=2021,
        help="Include papers published from this year (inclusive). Default: 2021.",
    )
    parser.add_argument(
        "--year-to",
        type=int,
        default=2026,
        help="Include papers published up to this year (inclusive). Default: 2026.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    library_id = extract_library_id(args.library)

    if not args.token:
        print(
            "Error: ADS API token is required. Set ADS_API_TOKEN or pass --token.",
            file=sys.stderr,
        )
        return 2
    if args.year_from > args.year_to:
        print("Error: --year-from must be <= --year-to.", file=sys.stderr)
        return 2

    client = AdsClient(token=args.token)
    try:
        metrics = compute_metrics(
            client,
            library_id,
            year_from=args.year_from,
            year_to=args.year_to,
        )
    except Exception as exc:  # noqa: BLE001
        print(f"Failed to compute metrics: {exc}", file=sys.stderr)
        return 1

    print(json.dumps(metrics, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
