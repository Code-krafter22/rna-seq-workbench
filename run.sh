#!/usr/bin/env bash
# Launch the RNA-seq DE web app.
set -e
cd "$(dirname "$0")"
exec streamlit run app/streamlit_app.py "$@"
