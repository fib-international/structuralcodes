selector_to_html = {"a[href=\"#installation\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">Installation<a class=\"headerlink\" href=\"#installation\" title=\"Link to this heading\">\u00b6</a></h1><p><code class=\"docutils literal notranslate\"><span class=\"pre\">structuralcodes</span></code> is compatible with Python 3.8, 3.9, 3.10, 3.11 and 3.12, and is installed by typing:</p>"}
skip_classes = ["headerlink", "sd-stretched-link"]

window.onload = function () {
    for (const [select, tip_html] of Object.entries(selector_to_html)) {
        const links = document.querySelectorAll(` ${select}`);
        for (const link of links) {
            if (skip_classes.some(c => link.classList.contains(c))) {
                continue;
            }

            tippy(link, {
                content: tip_html,
                allowHTML: true,
                arrow: true,
                placement: 'auto-start', maxWidth: 500, interactive: true,

            });
        };
    };
    console.log("tippy tips loaded!");
};
