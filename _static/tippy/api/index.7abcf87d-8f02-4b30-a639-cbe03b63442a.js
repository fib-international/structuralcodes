selector_to_html = {"a[href=\"materials/index.html#api-materials\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">Materials<a class=\"headerlink\" href=\"#materials\" title=\"Link to this heading\">\u00b6</a></h1>", "a[href=\"#api-reference\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">API reference<a class=\"headerlink\" href=\"#api-reference\" title=\"Link to this heading\">\u00b6</a></h1><p>This section contains descriptions of all modules, functions, classes and methods, as well as the structure of the library.</p>", "a[href=\"library_structure.html#api-structure\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">Library structure<a class=\"headerlink\" href=\"#library-structure\" title=\"Link to this heading\">\u00b6</a></h1><p>The content of StructuralCodes is organized in a hierarchy:</p>", "a[href=\"sections/index.html#api-sections\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">Sections<a class=\"headerlink\" href=\"#sections\" title=\"Link to this heading\">\u00b6</a></h1>", "a[href=\"#materials\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\"><a class=\"reference internal\" href=\"materials/index.html#api-materials\"><span class=\"std std-ref\">Materials</span></a><a class=\"headerlink\" href=\"#materials\" title=\"Link to this heading\">\u00b6</a></h2><p>Classes that represent materials from design codes.</p>", "a[href=\"geometry/index.html#api-geometry\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">Geometry<a class=\"headerlink\" href=\"#geometry\" title=\"Link to this heading\">\u00b6</a></h1>", "a[href=\"#sections\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\"><a class=\"reference internal\" href=\"sections/index.html#api-sections\"><span class=\"std std-ref\">Sections</span></a><a class=\"headerlink\" href=\"#sections\" title=\"Link to this heading\">\u00b6</a></h2><p>Classes that represent design sections.</p>", "a[href=\"codes/index.html#api-design-codes\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">Design codes<a class=\"headerlink\" href=\"#design-codes\" title=\"Link to this heading\">\u00b6</a></h1>", "a[href=\"#design-codes\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\"><a class=\"reference internal\" href=\"codes/index.html#api-design-codes\"><span class=\"std std-ref\">Design codes</span></a><a class=\"headerlink\" href=\"#design-codes\" title=\"Link to this heading\">\u00b6</a></h2><p>The design code equations are grouped in different modules.</p>", "a[href=\"#library-structure\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\"><a class=\"reference internal\" href=\"library_structure.html#api-structure\"><span class=\"std std-ref\">Library structure</span></a><a class=\"headerlink\" href=\"#library-structure\" title=\"Link to this heading\">\u00b6</a></h2><p>The content of the library is organized in a hierarchic manner.</p>", "a[href=\"#geometry\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\"><a class=\"reference internal\" href=\"geometry/index.html#api-geometry\"><span class=\"std std-ref\">Geometry</span></a><a class=\"headerlink\" href=\"#geometry\" title=\"Link to this heading\">\u00b6</a></h2><p>Classes for creating the geometry of design sections.</p>"}
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
