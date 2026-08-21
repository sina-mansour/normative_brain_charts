import re

def on_page_content(html, page, config, files):
    # Add target="_blank" to all external links
    html = re.sub(
        r'<a href="(https?://[^"]+)"',
        r'<a href="\1" target="_blank" rel="noopener"',
        html
    )
    return html
