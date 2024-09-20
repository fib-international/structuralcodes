selector_to_html = {"a[href=\"../quickstart/index.html#quickstart\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">Quickstart<a class=\"headerlink\" href=\"#quickstart\" title=\"Link to this heading\">\u00b6</a></h1><p>This example shows how to use <code class=\"docutils literal notranslate\"><span class=\"pre\">structuralcodes</span></code> to calculate the response of a rectangular reinforced concrete section. Follow the example step-by-step, or <a class=\"reference internal\" href=\"#quickstart-full-example\"><span class=\"std std-ref\">skip to the end</span></a> if you are in a hurry.</p><p>Import relevant functions and classes:</p>", "a[href=\"#examples\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">Examples<a class=\"headerlink\" href=\"#examples\" title=\"Link to this heading\">\u00b6</a></h1><p>We are working on a set of examples that demonstrate the use of the library. For now, have a look at the <a class=\"reference internal\" href=\"../quickstart/index.html#quickstart\"><span class=\"std std-ref\">quickstart</span></a>.</p>"}
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
