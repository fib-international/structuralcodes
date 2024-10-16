selector_to_html = {"a[href=\"#version-changes\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Version changes<a class=\"headerlink\" href=\"#version-changes\" title=\"Link to this heading\">\u00b6</a></h2><p>A <strong>patch</strong> version number increase, for example from <code class=\"docutils literal notranslate\"><span class=\"pre\">0.0.1</span></code> to <code class=\"docutils literal notranslate\"><span class=\"pre\">0.0.2</span></code>, occurs when:</p>", "a[href=\"#versioning\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">Versioning<a class=\"headerlink\" href=\"#versioning\" title=\"Link to this heading\">\u00b6</a></h1><h2>General<a class=\"headerlink\" href=\"#general\" title=\"Link to this heading\">\u00b6</a></h2><p>StructuralCodes does not yet have a stable API. This does not mean that you should not use the library for solving your design problems, but that there are more features to come, and that some of these new features might be implemented in ways that break existing behaviour.</p><p>Until then, StructuralCodes will use a custom versioning scheme that uses the <strong>minor</strong> version number for breaking changes, and the <strong>patch</strong> version number for other changes.</p>", "a[href=\"#general\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">General<a class=\"headerlink\" href=\"#general\" title=\"Link to this heading\">\u00b6</a></h2><p>StructuralCodes does not yet have a stable API. This does not mean that you should not use the library for solving your design problems, but that there are more features to come, and that some of these new features might be implemented in ways that break existing behaviour.</p><p>Until then, StructuralCodes will use a custom versioning scheme that uses the <strong>minor</strong> version number for breaking changes, and the <strong>patch</strong> version number for other changes.</p>"}
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
