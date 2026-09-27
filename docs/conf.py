# Sphinx configuration for penzene.readthedocs.io. Pages are Markdown (MyST); the Furo theme
# carries the app's white and teal colours.
project = "Penzene"
author = "James O'Brien"
copyright = "James O'Brien"

extensions = ["myst_parser", "sphinx_copybutton"]
source_suffix = {".md": "markdown"}
exclude_patterns = ["_build", "shortcuts.md"]
myst_heading_anchors = 3
myst_enable_extensions = ["colon_fence"]

html_theme = "furo"
html_title = "Penzene"
html_logo = "_static/logo.svg"
html_favicon = "_static/logo.svg"
html_static_path = ["_static"]
html_css_files = ["penzene.css"]
html_theme_options = {
    "source_repository": "https://github.com/JamesOBrien2/penzene/",
    "source_branch": "main",
    "source_directory": "docs/",
    "light_css_variables": {
        "color-brand-primary": "#0F6E56",
        "color-brand-content": "#0F6E56",
        "color-brand-visited": "#0F6E56",
        "color-background-primary": "#FFFFFF",
        "color-background-secondary": "#F3F5F6",
        "color-background-hover": "#E1F5EE",
        "color-background-border": "#DDE1E4",
        "color-foreground-primary": "#1D2125",
        "color-foreground-secondary": "#5E6770",
        "color-sidebar-background": "#F3F5F6",
        "color-sidebar-item-background--hover": "#E1F5EE",
        "color-code-background": "#FFFFFF",
        "color-admonition-background": "#FFFFFF",
        "color-admonition-title--note": "#0F6E56",
        "color-admonition-title-background--note": "#E1F5EE",
    },
    "dark_css_variables": {
        "color-brand-primary": "#4CC9A0",
        "color-brand-content": "#4CC9A0",
        "color-brand-visited": "#4CC9A0",
        "color-background-primary": "#15171A",
        "color-background-secondary": "#1C1F23",
        "color-background-hover": "#0B3B30",
        "color-background-border": "#353A40",
        "color-foreground-primary": "#E8EAED",
        "color-foreground-secondary": "#A3AAB2",
        "color-sidebar-background": "#1C1F23",
        "color-sidebar-item-background--hover": "#0B3B30",
        "color-code-background": "#24282D",
        "color-admonition-background": "#24282D",
        "color-admonition-title--note": "#4CC9A0",
        "color-admonition-title-background--note": "#0B3B30",
    },
}
copybutton_prompt_text = r"\$ |>>> "
copybutton_prompt_is_regexp = True

# ```{feature} icon-name
# :title: Formula and mass
# Markdown body.
# ```
# A What's New style card: a teal icon badge (a Tabler SVG from _static/icons, inlined so it
# takes the theme's colour), a title and the text.
from pathlib import Path
from docutils import nodes
from docutils.parsers.rst import Directive, directives

ICONS = Path(__file__).parent / "_static" / "icons"


class Feature(Directive):
    required_arguments = 1
    option_spec = {"title": directives.unchanged_required}
    has_content = True

    def run(self):
        svg = (ICONS / f"{self.arguments[0]}.svg").read_text(encoding="utf-8")
        card = nodes.container(classes=["feature"])
        card += nodes.raw("", f'<span class="badge">{svg}</span>', format="html")
        body = nodes.container(classes=["feature-body"])
        body += nodes.paragraph("", "", nodes.strong(text=self.options["title"]), classes=["feature-title"])
        self.state.nested_parse(self.content, self.content_offset, body)
        card += body
        return [card]


def setup(app):
    app.add_directive("feature", Feature)
