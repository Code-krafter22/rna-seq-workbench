"""Presentation layer for the RNA-seq workbench."""

from __future__ import annotations

import html

import streamlit as st

FONT_URL = (
    "https://fonts.googleapis.com/css2?"
    "family=IBM+Plex+Mono:wght@400;500&"
    "family=IBM+Plex+Sans:wght@400;500;600&"
    "family=Newsreader:opsz,wght@6..72,400;6..72,500&display=swap"
)

CSS = f"""
<style>
@import url('{FONT_URL}');

:root {{
  --paper: #FAFBF8;
  --surface: #FFFFFF;
  --surface-2: #F2F5EF;
  --ink: #14211C;
  --ink-soft: #33443B;
  --muted: #5C6C64;
  --muted-2: #7C8C84;
  --line: #DEE5DC;
  --line-strong: #CBD6CE;
  --accent: #0E6B54;
  --accent-strong: #0A5743;
  --accent-soft: #E9F3ED;
  --mint: #8FD9BC;
  --amber: #9A6414;
  --red: #A33A2E;
  --sans: "IBM Plex Sans", ui-sans-serif, -apple-system, "Segoe UI", sans-serif;
  --serif: "Newsreader", ui-serif, "Iowan Old Style", Georgia, serif;
  --mono: "IBM Plex Mono", ui-monospace, SFMono-Regular, Menlo, monospace;
  --sidebar: #0F1C17;
  --sidebar-line: #2A3D34;
  --sidebar-muted: #9DB0A5;
}}

html, body, .stApp, [data-testid="stAppViewContainer"] {{
  background: var(--paper) !important;
  font-family: var(--sans);
}}

html {{ -webkit-font-smoothing: antialiased; }}

::selection {{ background: #CDE8DB; color: var(--ink); }}

#MainMenu, footer, [data-testid="stAppDeployButton"] {{ display: none !important; }}

[data-testid="stHeader"] {{ background: transparent; }}
[data-testid="stToolbar"] {{ right: 0.75rem; top: 0.35rem; }}

[data-testid="stMainBlockContainer"] {{
  max-width: 1200px;
  padding: 1.5rem 2.5rem 3rem;
}}

.stApp p, .stApp li, .stApp label, .stApp input, .stApp textarea,
.stApp button, .stApp select, .stApp [data-testid="stMarkdownContainer"] {{
  font-family: var(--sans);
}}

.stApp h1, .stApp h2 {{
  font-family: var(--serif);
  color: var(--ink);
  letter-spacing: -0.02em;
  font-weight: 500;
}}

.stApp h3, .stApp h4, .stApp h5 {{
  font-family: var(--sans);
  color: var(--ink);
  font-weight: 600;
}}

.stApp h3 {{
  font-size: 1.05rem;
  letter-spacing: -0.005em;
  margin: 1.35rem 0 0.2rem;
}}

.stApp a {{ color: var(--accent-strong); }}
.stApp a:hover {{ color: var(--accent); }}

[data-testid="stCaptionContainer"], [data-testid="stCaptionContainer"] p {{
  color: var(--muted);
  font-size: 0.8rem;
  line-height: 1.55;
}}

[data-testid="stMarkdownContainer"] p {{ line-height: 1.62; }}

hr, [data-testid="stDivider"] hr {{
  border-color: var(--line);
  margin: 1.5rem 0;
}}

/* sidebar */
[data-testid="stSidebar"] {{
  background: var(--sidebar) !important;
  border-right: 1px solid #21342B;
}}

[data-testid="stSidebarContent"] {{ background: var(--sidebar); }}

[data-testid="stSidebarUserContent"] {{
  padding: 0.5rem 1.15rem 1.5rem;
}}

[data-testid="stSidebarHeader"] {{
  padding: 0.6rem 0.9rem 0;
}}

[data-testid="stSidebarCollapseButton"] button,
[data-testid="stSidebarCollapseButton"] span {{
  color: #C9D6CD !important;
}}

.od-brand {{
  display: flex;
  align-items: center;
  gap: 0.7rem;
  padding: 0.35rem 0.15rem 1.35rem;
}}

.od-brand-mark {{
  width: 2.4rem;
  height: 2.4rem;
  flex: none;
  display: grid;
  place-items: center;
  border-radius: 0.6rem;
  color: var(--mint);
  background: rgba(143, 217, 188, 0.1);
  border: 1px solid rgba(143, 217, 188, 0.22);
}}

.od-brand-mark svg {{ width: 1.4rem; height: 1.4rem; }}

.od-brand-name {{
  margin: 0;
  color: #F3F8F4;
  font-family: var(--serif);
  font-size: 1.08rem;
  line-height: 1.15;
  letter-spacing: -0.01em;
}}

.od-brand-sub {{
  margin: 0.2rem 0 0;
  color: var(--sidebar-muted);
  font-family: var(--mono);
  font-size: 0.55rem;
  letter-spacing: 0.14em;
  text-transform: uppercase;
  line-height: 1.4;
}}

.od-nav-label {{
  margin: 0 0 0.4rem 0.8rem;
  color: #7E958A;
  font-family: var(--mono);
  font-size: 0.56rem;
  letter-spacing: 0.16em;
  text-transform: uppercase;
}}

[data-testid="stSidebar"] [role="radiogroup"] {{
  gap: 0.12rem;
  counter-reset: odstep;
}}

[data-testid="stSidebar"] [role="radiogroup"] > label {{
  counter-increment: odstep;
  display: flex;
  align-items: center;
  gap: 0.55rem;
  min-height: 2.75rem;
  padding: 0.55rem 0.7rem;
  border-radius: 0.5rem;
  color: #C9D6CD !important;
  cursor: pointer;
  transition: background-color 160ms ease, color 160ms ease;
}}

[data-testid="stSidebar"] [role="radiogroup"] > label::before {{
  content: counter(odstep, decimal-leading-zero);
  flex: none;
  color: #8FA697;
  font-family: var(--mono);
  font-size: 0.6rem;
  letter-spacing: 0.06em;
  transition: color 160ms ease;
}}

[data-testid="stSidebar"] [role="radiogroup"] > label p {{
  color: inherit !important;
  font-size: 0.84rem;
  font-weight: 500;
  line-height: 1.3;
}}

[data-testid="stSidebar"] [role="radiogroup"] > label:hover {{
  background: rgba(255, 255, 255, 0.055);
}}

[data-testid="stSidebar"] [role="radiogroup"] > label:has(input:checked) {{
  background: rgba(143, 217, 188, 0.12);
  box-shadow: inset 2px 0 0 var(--mint);
  color: #F7FBF8 !important;
}}

[data-testid="stSidebar"] [role="radiogroup"] > label:has(input:checked)::before {{
  color: var(--mint);
}}

[data-testid="stSidebar"] [role="radiogroup"] > label:has(input:focus-visible) {{
  outline: 2px solid var(--mint);
  outline-offset: -2px;
}}

.od-status {{
  margin-top: 1.35rem;
  padding: 0.85rem 0.85rem 0.55rem;
  border: 1px solid var(--sidebar-line);
  border-radius: 0.75rem;
  background: rgba(255, 255, 255, 0.035);
}}

.od-status-title {{
  margin: 0 0 0.45rem;
  color: #7E958A;
  font-family: var(--mono);
  font-size: 0.56rem;
  letter-spacing: 0.16em;
  text-transform: uppercase;
}}

.od-status-row {{
  display: flex;
  align-items: baseline;
  justify-content: space-between;
  gap: 0.7rem;
  padding: 0.42rem 0;
  border-top: 1px solid rgba(255, 255, 255, 0.06);
}}

.od-status-row:first-of-type {{ border-top: 0; }}

.od-status-k {{
  color: var(--sidebar-muted);
  font-size: 0.73rem;
  line-height: 1.4;
}}

.od-status-v {{
  max-width: 10.5rem;
  color: #E8EFE9;
  font-size: 0.73rem;
  font-weight: 500;
  line-height: 1.4;
  overflow-wrap: anywhere;
  text-align: right;
}}

.od-status-v.is-ok {{ color: var(--mint); }}
.od-status-v.is-idle {{ color: #93A89D; font-weight: 400; }}

.od-sidefoot {{
  margin-top: 0.85rem;
  color: #8CA196;
  font-size: 0.65rem;
  line-height: 1.55;
}}

/* institution band (USF colors: green #006747, gold #CFC493) */
.od-inst {{
  display: flex;
  align-items: center;
  gap: 0.75rem;
  margin: 0 0 0.9rem;
  padding: 0.6rem 0.9rem;
  border-left: 4px solid #CFC493;
  border-radius: 0.6rem;
  background: #006747;
  color: #FFFFFF;
}}
.od-inst-mark {{
  display: inline-grid;
  place-items: center;
  min-width: 2.6rem;
  height: 2.6rem;
  border-radius: 0.45rem;
  background: #CFC493;
  color: #00402C;
  font-family: var(--serif);
  font-weight: 600;
  font-size: 1.05rem;
  letter-spacing: 0.04em;
}}
.od-inst-text {{ display: flex; flex-direction: column; line-height: 1.2; }}
.od-inst-text b {{ font-weight: 600; font-size: 1rem; }}
.od-inst-text span {{ font-size: 0.85rem; color: #DCEFE5; letter-spacing: 0.02em; }}

/* hero */
.od-hero {{
  display: grid;
  grid-template-columns: minmax(0, 1.5fr) minmax(0, 0.85fr);
  gap: 2.5rem;
  align-items: end;
  margin: 0.25rem 0 0.5rem;
  padding: 2.3rem 2.6rem 2.1rem;
  border: 1px solid var(--line);
  border-radius: 1.1rem;
  background: var(--surface);
  box-shadow: 0 1px 2px rgba(16, 32, 26, 0.035);
}}

.od-eyebrow {{
  display: flex;
  align-items: center;
  gap: 0.55rem;
  margin: 0;
  color: var(--accent);
  font-family: var(--mono);
  font-size: 0.62rem;
  font-weight: 500;
  letter-spacing: 0.15em;
  line-height: 1.4;
  text-transform: uppercase;
}}

.od-eyebrow::before {{
  content: "";
  width: 1.1rem;
  height: 1px;
  background: var(--accent);
}}

.od-hero h1 {{
  margin: 0.7rem 0 0.75rem;
  font-size: clamp(2.35rem, 4vw, 3.5rem);
  line-height: 1.02;
  letter-spacing: -0.028em;
}}

.od-lede {{
  max-width: 54ch;
  margin: 0;
  color: var(--muted);
  font-size: 0.96rem;
  line-height: 1.65;
}}

.od-hero-foot {{
  display: flex;
  flex-wrap: wrap;
  align-items: center;
  justify-content: space-between;
  gap: 0.75rem 1.25rem;
  margin-top: 1.5rem;
}}

.od-chips {{ display: flex; flex-wrap: wrap; gap: 0.35rem; }}

.od-chip {{
  padding: 0.22rem 0.5rem;
  border: 1px solid var(--line);
  border-radius: 999px;
  background: var(--surface-2);
  color: var(--ink-soft);
  font-family: var(--mono);
  font-size: 0.62rem;
  letter-spacing: 0.02em;
  line-height: 1.5;
}}

.od-hero-status {{
  display: inline-flex;
  align-items: center;
  gap: 0.45rem;
  color: var(--muted);
  font-size: 0.76rem;
  line-height: 1.4;
}}

.od-dot {{
  width: 0.42rem;
  height: 0.42rem;
  flex: none;
  border-radius: 50%;
  background: var(--line-strong);
}}

.od-dot.is-live {{ background: var(--accent); box-shadow: 0 0 0 3px var(--accent-soft); }}

.od-hero-track {{
  align-self: stretch;
  display: flex;
  align-items: flex-end;
  min-width: 0;
  padding-top: 0.5rem;
}}

.od-track-bars {{
  display: flex;
  align-items: flex-end;
  gap: 3px;
  width: 100%;
  height: 9.4rem;
  padding-bottom: 1px;
  border-bottom: 1px solid var(--line-strong);
}}

.od-track-bars i {{
  display: block;
  flex: 1 1 0;
  min-width: 2px;
  border-radius: 2px 2px 0 0;
  background: #DFE9E1;
}}

.od-track-bars i.is-mid {{ background: var(--mint); }}
.od-track-bars i.is-hot {{ background: var(--accent); }}

/* step header */
.od-step-head {{
  margin: 1.9rem 0 1.25rem;
  padding-bottom: 0.95rem;
  border-bottom: 1px solid var(--line);
}}

.od-step-kicker {{
  display: flex;
  gap: 0.4rem;
  margin: 0;
  color: var(--accent);
  font-family: var(--mono);
  font-size: 0.62rem;
  font-weight: 500;
  letter-spacing: 0.14em;
  line-height: 1.4;
  text-transform: uppercase;
}}

.od-step-kicker span {{ color: var(--muted); }}

.od-step-head h2 {{
  margin: 0.3rem 0 0.35rem;
  font-size: clamp(1.9rem, 2.6vw, 2.45rem);
  line-height: 1.08;
}}

.od-step-lede {{
  max-width: 74ch;
  margin: 0;
  color: var(--muted);
  font-size: 0.93rem;
  line-height: 1.6;
}}

/* navigation */
.od-nav-meta {{
  display: flex;
  align-items: center;
  gap: 0.9rem;
  margin: 2rem 0 0.65rem;
}}

.od-nav-step {{
  color: var(--muted);
  font-family: var(--mono);
  font-size: 0.6rem;
  letter-spacing: 0.1em;
  line-height: 1.4;
  text-transform: uppercase;
  white-space: nowrap;
}}

.st-key-nav_row [data-testid="stHorizontalBlock"] {{
  max-width: 34rem;
}}

.od-nav-track {{
  flex: 1;
  height: 2px;
  overflow: hidden;
  border-radius: 2px;
  background: var(--line);
}}

.od-nav-track i {{
  display: block;
  height: 100%;
  border-radius: 2px;
  background: var(--accent);
}}

.od-footer {{
  display: flex;
  flex-wrap: wrap;
  justify-content: space-between;
  gap: 0.4rem 1.5rem;
  margin-top: 1.9rem;
  padding-top: 0.9rem;
  border-top: 1px solid var(--line);
  color: var(--muted);
  font-size: 0.68rem;
  line-height: 1.55;
}}

/* auth */
.od-auth {{
  max-width: 30rem;
  margin: 10vh auto 1.2rem;
  padding: 2rem 2.1rem;
  border: 1px solid var(--line);
  border-radius: 1rem;
  background: var(--surface);
}}

.od-auth h2 {{
  margin: 0.5rem 0 0.4rem;
  font-size: 2rem;
  line-height: 1.1;
}}

.od-auth .od-lede {{ font-size: 0.9rem; }}

/* native widgets */
[data-testid="stButton"] > button, [data-testid="stDownloadButton"] > button,
.stButton > button, .stDownloadButton > button {{
  min-height: 2.75rem;
  border-radius: 8px;
  font-size: 0.84rem;
  font-weight: 500;
  letter-spacing: 0.005em;
  transition: background-color 160ms ease, border-color 160ms ease,
              color 160ms ease, transform 160ms ease;
}}

[data-testid="stButton"] > button:active,
[data-testid="stDownloadButton"] > button:active,
.stButton > button:active, .stDownloadButton > button:active {{
  transform: translateY(1px);
}}

[data-testid="stButton"] > button:focus-visible,
[data-testid="stDownloadButton"] > button:focus-visible,
.stButton > button:focus-visible, .stDownloadButton > button:focus-visible,
.stApp a:focus-visible {{
  outline: 3px solid rgba(14, 107, 84, 0.28);
  outline-offset: 2px;
}}

button[kind="primary"], [data-testid="stBaseButton-primary"] {{
  border: 1px solid var(--accent) !important;
  background: var(--accent) !important;
  color: #FFFFFF !important;
}}

button[kind="primary"]:hover, [data-testid="stBaseButton-primary"]:hover {{
  border-color: var(--accent-strong) !important;
  background: var(--accent-strong) !important;
  color: #FFFFFF !important;
}}

button[kind="secondary"], [data-testid="stBaseButton-secondary"] {{
  border: 1px solid var(--line-strong) !important;
  background: var(--surface) !important;
  color: var(--ink-soft) !important;
}}

button[kind="secondary"]:hover, [data-testid="stBaseButton-secondary"]:hover {{
  border-color: var(--accent) !important;
  background: var(--accent-soft) !important;
  color: var(--accent-strong) !important;
}}

button[kind="tertiary"], [data-testid="stBaseButton-tertiary"] {{
  color: var(--accent-strong) !important;
}}

[data-testid="stWidgetLabel"] p {{
  color: var(--ink-soft);
  font-size: 0.82rem;
  font-weight: 600;
  line-height: 1.4;
}}

[data-testid="stMain"] [role="radiogroup"] > label,
[data-testid="stMain"] [data-testid="stCheckbox"] label {{
  font-size: 0.85rem;
  color: var(--ink-soft);
}}

[data-testid="stMain"] [role="radiogroup"] > label > div:first-of-type {{
  border-color: var(--line-strong) !important;
}}

[data-testid="stMain"] [role="radiogroup"] > label:has(input:checked) > div:first-of-type {{
  border-color: var(--accent) !important;
}}

[data-testid="stMain"] [role="radiogroup"] > label:has(input:checked) > div:first-of-type > div {{
  background-color: var(--accent) !important;
}}

[data-testid="stSidebar"] [role="radiogroup"] > label > div:first-of-type {{
  border-color: #51675C !important;
}}

[data-testid="stSidebar"] [role="radiogroup"] > label:has(input:checked) > div:first-of-type {{
  border-color: var(--mint) !important;
}}

[data-testid="stSidebar"] [role="radiogroup"] > label:has(input:checked) > div:first-of-type > div {{
  background-color: var(--mint) !important;
}}

[data-testid="stCheckbox"] label:has(input:checked) > div:first-of-type {{
  border-color: var(--accent) !important;
  background-color: var(--accent) !important;
}}

[data-baseweb="input"]:focus-within,
[data-baseweb="textarea"]:focus-within,
[data-baseweb="select"]:focus-within {{
  border-color: var(--accent) !important;
  box-shadow: 0 0 0 3px rgba(14, 107, 84, 0.12);
}}

[data-testid="stSlider"] [role="slider"] {{
  border-color: var(--accent) !important;
  background-color: var(--accent) !important;
}}

[data-testid="stFileUploaderDropzone"] {{
  padding: 1rem;
  border: 1px dashed var(--line-strong);
  border-radius: 10px;
  background: #FBFCF9;
  transition: border-color 160ms ease, background-color 160ms ease;
}}

[data-testid="stFileUploaderDropzone"]:hover {{
  border-color: var(--accent);
  background: var(--accent-soft);
}}

[data-testid="stFileUploaderDropzoneInstructions"] span,
[data-testid="stFileUploaderDropzoneInstructions"] small {{
  color: var(--muted);
  font-size: 0.78rem;
}}

[data-testid="stFileUploaderDropzone"] button {{
  border-radius: 8px;
  font-size: 0.8rem;
}}

[data-testid="stAlertContainer"] {{
  border: 1px solid rgba(20, 33, 28, 0.06);
  border-radius: 10px;
  padding: 0.7rem 0.85rem;
}}

[data-testid="stAlertContainer"] p {{
  font-size: 0.84rem;
  line-height: 1.55;
}}

[data-testid="stMetric"] {{
  padding: 0.9rem 1rem;
  border: 1px solid var(--line);
  border-radius: 10px;
  background: var(--surface);
}}

[data-testid="stMetricLabel"] p {{
  color: var(--muted);
  font-size: 0.64rem;
  font-weight: 600;
  letter-spacing: 0.09em;
  line-height: 1.4;
  text-transform: uppercase;
}}

[data-testid="stMetricValue"] {{
  color: var(--ink);
  font-family: var(--mono);
  font-size: 1.65rem;
  font-weight: 500;
  letter-spacing: -0.02em;
  line-height: 1.25;
}}

[data-testid="stMetricDelta"] {{
  color: var(--muted);
  font-family: var(--mono);
  font-size: 0.66rem;
}}

[data-testid="stTabs"] [data-baseweb="tab-list"] {{
  gap: 0.3rem;
  border-bottom: 1px solid var(--line);
}}

[data-testid="stTabs"] button[data-baseweb="tab"] {{
  padding: 0.45rem 0.7rem;
  color: var(--muted);
  font-size: 0.82rem;
}}

[data-testid="stTabs"] button[aria-selected="true"] {{ color: var(--accent); }}

[data-testid="stTabs"] [data-baseweb="tab-highlight"] {{ background: var(--accent); }}

[data-testid="stExpander"] details {{
  border: 1px solid var(--line);
  border-radius: 10px;
  background: var(--surface);
  overflow: hidden;
}}

[data-testid="stExpander"] summary {{
  padding: 0.65rem 0.85rem;
  color: var(--ink-soft);
  font-size: 0.82rem;
  font-weight: 600;
}}

[data-testid="stExpander"] summary:hover {{ color: var(--accent); }}

[data-testid="stExpander"] [data-testid="stExpanderDetails"] {{
  padding: 0.1rem 0.9rem 0.9rem;
}}

[data-testid="stDataFrame"] {{
  overflow: hidden;
  border: 1px solid var(--line);
  border-radius: 10px;
}}

[data-testid="stImage"] img {{
  display: block;
  width: 100%;
  height: auto;
  border: 1px solid var(--line);
  border-radius: 8px;
  background: var(--surface);
}}

[data-testid="stImage"] [data-testid="stImageCaption"] {{
  color: var(--muted);
  font-family: var(--mono);
  font-size: 0.66rem;
}}

[data-testid="stCodeBlock"] code, [data-testid="stCode"] code {{
  font-family: var(--mono);
  font-size: 0.76rem;
}}

[data-testid="stPlotlyChart"], [data-testid="stArrowVegaLiteChart"] {{
  border: 1px solid var(--line);
  border-radius: 10px;
  background: var(--surface);
}}

@media (max-width: 1080px) {{
  .od-hero {{
    grid-template-columns: minmax(0, 1fr);
    gap: 1.4rem;
    padding: 2rem 1.9rem 1.7rem;
  }}
  .od-hero-track {{ display: none; }}
}}

@media (max-width: 760px) {{
  [data-testid="stMainBlockContainer"] {{
    padding: 1rem 1rem 2.5rem;
  }}
  .od-hero {{
    padding: 1.5rem 1.25rem 1.35rem;
    border-radius: 0.9rem;
  }}
  .od-hero h1 {{ font-size: 2.25rem; }}
  .od-lede {{ font-size: 0.92rem; }}
  .od-step-head {{ margin-top: 1.4rem; }}
  .od-step-head h2 {{ font-size: 1.8rem; }}
  .od-footer {{ flex-direction: column; }}
}}

@media (prefers-reduced-motion: reduce) {{
  *, *::before, *::after {{
    animation-duration: 0.01ms !important;
    animation-iteration-count: 1 !important;
    transition-duration: 0.01ms !important;
    scroll-behavior: auto !important;
  }}
}}
</style>
"""

_BRAND_MARK = """
<svg viewBox="0 0 32 32" fill="none" aria-hidden="true">
  <path d="M10 4c12 5 12 19 0 24" stroke="currentColor" stroke-width="1.8"
        stroke-linecap="round"/>
  <path d="M22 4c-12 5-12 19 0 24" stroke="currentColor" stroke-width="1.8"
        stroke-linecap="round"/>
  <path d="M11.5 11h9M11.5 21h9" stroke="currentColor" stroke-width="1.6"
        stroke-linecap="round"/>
</svg>
"""

_HERO_BARS = [
    22, 34, 18, 46, 30, 58, 40, 26, 52, 70, 44, 32, 62, 48, 24, 38,
    56, 76, 42, 30, 50, 36, 64, 46, 28, 54, 40, 66, 34, 48, 58, 42,
]


def inject_css() -> None:
    st.markdown(CSS, unsafe_allow_html=True)


def _esc(value: object) -> str:
    return html.escape(str(value), quote=True)


def contrast_label(contrast: str | None) -> str | None:
    if not contrast:
        return None
    parts = [part.strip() for part in contrast.split(",") if part.strip()]
    if len(parts) == 2:
        return f"{parts[0]} vs {parts[1]}"
    return contrast.strip()


def hero(
    sample_count: int | None,
    counts_loaded: bool,
    run_ready: bool,
    contrast: str | None,
) -> None:
    state_bits: list[str] = []
    if sample_count is not None:
        state_bits.append(f"{sample_count} samples")
    label = contrast_label(contrast)
    if label:
        state_bits.append(label)
    if run_ready:
        state_bits.append("results ready")
    elif counts_loaded:
        state_bits.append("counts loaded")
    state = " · ".join(state_bits) if state_bits else "No dataset loaded"
    live = " is-live" if (counts_loaded or sample_count or run_ready) else ""

    hot = {9, 17, 22, 27}
    mid = {5, 12, 18, 26}
    bars = "".join(
        '<i class="{cls}" style="height:{height}%"></i>'.format(
            cls="is-hot" if index in hot else ("is-mid" if index in mid else ""),
            height=height,
        )
        for index, height in enumerate(_HERO_BARS)
    )

    st.markdown(
        f"""
<section class="od-hero">
  <div class="od-hero-body">
    <p class="od-eyebrow">Bulk RNA-seq · differential expression</p>
    <h1>From raw counts to a defensible gene list.</h1>
    <p class="od-lede">Run DESeq2, edgeR and limma-voom side by side on a GEO
    matrix or your own counts, then check batch, read the results, and follow
    a gene into the literature.</p>
    <div class="od-hero-foot">
      <div class="od-chips">
        <span class="od-chip">DESeq2</span>
        <span class="od-chip">edgeR</span>
        <span class="od-chip">limma-voom</span>
      </div>
      <span class="od-hero-status"><span class="od-dot{live}"></span>{_esc(state)}</span>
    </div>
  </div>
  <div class="od-hero-track" aria-hidden="true">
    <div class="od-track-bars">{bars}</div>
  </div>
</section>
""",
        unsafe_allow_html=True,
    )


def step_header(number: int, total: int, title: str, lede: str) -> None:
    st.markdown(
        f"""
<div class="od-step-head">
  <p class="od-step-kicker">Step {number:02d}<span>/ {total:02d}</span></p>
  <h2>{_esc(title)}</h2>
  <p class="od-step-lede">{_esc(lede)}</p>
</div>
""",
        unsafe_allow_html=True,
    )


def sidebar_brand() -> None:
    st.sidebar.markdown(
        f"""
<div class="od-brand">
  <div class="od-brand-mark">{_BRAND_MARK}</div>
  <div>
    <p class="od-brand-name">RNA-seq Workbench</p>
    <p class="od-brand-sub">Differential expression</p>
  </div>
</div>
<p class="od-nav-label">Workflow</p>
""",
        unsafe_allow_html=True,
    )


def _status_row(label: str, value: str, ok: bool, idle: bool = False) -> str:
    cls = "is-ok" if ok else ("is-idle" if idle else "")
    return (
        '<div class="od-status-row">'
        f'<span class="od-status-k">{_esc(label)}</span>'
        f'<span class="od-status-v {cls}" title="{_esc(value)}">{_esc(value)}</span>'
        "</div>"
    )


def sidebar_status(
    counts_label: str,
    samples_label: str,
    results_label: str,
    contrast_label: str | None,
    batch_label: str | None = None,
) -> None:
    rows = [
        _status_row("Counts matrix", counts_label, counts_label == "Loaded",
                    idle=counts_label != "Loaded"),
        _status_row("Sample sheet", samples_label, samples_label != "Not loaded",
                    idle=samples_label == "Not loaded"),
        _status_row("Results", results_label, results_label == "Ready",
                    idle=results_label != "Ready"),
        _status_row("Contrast", contrast_label or "Not set",
                    bool(contrast_label), idle=not contrast_label),
    ]
    if batch_label:
        rows.append(_status_row("Batch", batch_label, True))

    st.sidebar.markdown(
        """
<div class="od-status">
  <p class="od-status-title">Workspace</p>
"""
        + "".join(rows)
        + """
</div>
<div class="od-sidefoot">Local Streamlit app · R/Bioconductor engine</div>
""",
        unsafe_allow_html=True,
    )


def nav_progress(index: int, total: int, label: str) -> None:
    percent = round((index + 1) / total * 100, 1)
    st.markdown(
        f"""
<div class="od-nav-meta">
  <span class="od-nav-step">Step {index + 1:02d} of {total:02d} · {_esc(label)}</span>
  <span class="od-nav-track"><i style="width:{percent}%"></i></span>
</div>
""",
        unsafe_allow_html=True,
    )


def footer() -> None:
    st.markdown(
        """
<div class="od-footer">
  <span>Bulk RNA-seq workbench · DESeq2 / edgeR / limma-voom</span>
  <span>GEO fetch, KEGG/Reactome enrichment and literature search use the network.</span>
</div>
""",
        unsafe_allow_html=True,
    )


def auth_header() -> None:
    st.markdown(
        """
<section class="od-auth">
  <p class="od-eyebrow">Restricted workspace</p>
  <h2>Sign in to the workbench</h2>
  <p class="od-lede">This instance is protected by a shared password. Enter it
  to load the RNA-seq workflow.</p>
</section>
""",
        unsafe_allow_html=True,
    )


def institution_band() -> None:
    st.markdown(
        """
<div class="od-inst" role="banner">
  <span class="od-inst-mark" aria-hidden="true">USF</span>
  <span class="od-inst-text">
    <b>University of South Florida</b>
    <span>College of Nursing</span>
  </span>
</div>
""",
        unsafe_allow_html=True,
    )


def scroll_to_top(step: str) -> None:
    """Jump to the top after a step change (Streamlit keeps scroll position across
    reruns). `step` is in the markup so each change mounts a fresh script."""
    st.html(
        f"""<script data-step="{_esc(step)}">
setTimeout(() => {{
  document.querySelectorAll('[data-testid="stMain"], [data-testid="stAppViewContainer"], section.main')
    .forEach((el) => el.scrollTo({{ top: 0 }}));
  window.scrollTo(0, 0);
}}, 60);
</script>""",
        unsafe_allow_javascript=True,
    )
