selector_to_html = {"a[href=\"#theory-reference\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">Theory reference<a class=\"headerlink\" href=\"#theory-reference\" title=\"Link to this heading\">\u00b6</a></h1><p>We are working on a description of the theoretical basis for the implementations in the library. Check back regularly for updates.</p>"}
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
