#!/usr/bin/env python3
"""Dump the Python scenarios as JSON, for tools/check_scenarios.js.

    python3 tools/dump_scenarios.py /tmp/scen.json
"""

from __future__ import annotations

import dataclasses
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from ventsim.scenarios import SCENARIOS  # noqa: E402


def main() -> int:
    out = {}
    for s in SCENARIOS:
        out[s.key] = {
            "title": s.title,
            "adapt_hours": s.adapt_hours,
            "story": s.story,
            "try": s.try_,
            "expect": s.expect,
            "patient": {k: v for k, v in dataclasses.asdict(s.patient).items()
                        if k != "lung"},
            "vent": dataclasses.asdict(s.vent),
        }
    dest = Path(sys.argv[1] if len(sys.argv) > 1 else "/tmp/scen.json")
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text(json.dumps(out, indent=1) + "\n")
    print(f"{len(out)} scenarios -> {dest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
