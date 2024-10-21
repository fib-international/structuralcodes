selector_to_html = {"a[href=\"#structuralcodes.codes.get_design_codes\"]": "<dt class=\"sig sig-object py\" id=\"structuralcodes.codes.get_design_codes\">\n<span class=\"sig-prename descclassname\"><span class=\"pre\">structuralcodes.codes.</span></span><span class=\"sig-name descname\"><span class=\"pre\">get_design_codes</span></span><span class=\"sig-paren\">(</span><span class=\"sig-paren\">)</span> <span class=\"sig-return\"><span class=\"sig-return-icon\">\u2192</span> <span class=\"sig-return-typehint\"><span class=\"pre\">List</span><span class=\"p\"><span class=\"pre\">[</span></span><span class=\"pre\">str</span><span class=\"p\"><span class=\"pre\">]</span></span></span></span><a class=\"reference external\" href=\"https://github.com/fib-international/structuralcodes/tree/main/structuralcodes\\codes\\__init__.py#L66-L68\"><span class=\"viewcode-link\"><span class=\"pre\">[source]</span></span></a></dt><dd><p>Get a list of the available design codes.</p></dd>", "a[href=\"#structuralcodes.codes.set_design_code\"]": "<dt class=\"sig sig-object py\" id=\"structuralcodes.codes.set_design_code\">\n<span class=\"sig-prename descclassname\"><span class=\"pre\">structuralcodes.codes.</span></span><span class=\"sig-name descname\"><span class=\"pre\">set_design_code</span></span><span class=\"sig-paren\">(</span><em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">design_code</span></span><span class=\"p\"><span class=\"pre\">:</span></span><span class=\"w\"> </span><span class=\"n\"><span class=\"pre\">str</span><span class=\"w\"> </span><span class=\"p\"><span class=\"pre\">|</span></span><span class=\"w\"> </span><span class=\"pre\">ModuleType</span><span class=\"w\"> </span><span class=\"p\"><span class=\"pre\">|</span></span><span class=\"w\"> </span><span class=\"pre\">None</span></span><span class=\"w\"> </span><span class=\"o\"><span class=\"pre\">=</span></span><span class=\"w\"> </span><span class=\"default_value\"><span class=\"pre\">None</span></span></em><span class=\"sig-paren\">)</span> <span class=\"sig-return\"><span class=\"sig-return-icon\">\u2192</span> <span class=\"sig-return-typehint\"><span class=\"pre\">None</span></span></span><a class=\"reference external\" href=\"https://github.com/fib-international/structuralcodes/tree/main/structuralcodes\\codes\\__init__.py#L33-L63\"><span class=\"viewcode-link\"><span class=\"pre\">[source]</span></span></a></dt><dd><p>Set the current design code globally.</p></dd>", "a[href=\"#general\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">General<a class=\"headerlink\" href=\"#general\" title=\"Link to this heading\">\u00b6</a></h1>"}
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
