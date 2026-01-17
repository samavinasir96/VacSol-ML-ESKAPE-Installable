import os
os.environ.setdefault(
    "DJANGO_SETTINGS_MODULE",
    "vacsol_ml.config.settings"
)

import argparse
import sys
from pathlib import Path
from django.core.management import execute_from_command_line



def main():
    parser = argparse.ArgumentParser(
        prog="vacsol-ml",
        description="VacSol-ML(ESKAPE) : ML-based vaccine target discovery targetting ESKAPE pathogens"
    )

    sub = parser.add_subparsers(dest="cmd", required=True)

    sub.add_parser("check", help="Check VacSol-ML installation")
    sub.add_parser("web", help="Run Django development server")
    sub.add_parser("migrate")

    args = parser.parse_args()

    if args.cmd == "check":
        import vacsol_ml
        print("✅ VacSol-ML installed")
        print("📦 Location:", vacsol_ml.__file__)

    elif args.cmd == "web":
        execute_from_command_line(
            ["manage.py", "runserver"]
        )

    elif args.cmd == "migrate":
        execute_from_command_line(
            ["manage.py", "migrate"]
        )

    else:
        parser.print_help()


if __name__ == "__main__":
    main()
