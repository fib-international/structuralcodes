selector_to_html = {"a[href=\"#structuralcodes.codes.mc2010.epsud\"]": "<dt class=\"sig sig-object py\" id=\"structuralcodes.codes.mc2010.epsud\">\n<span class=\"sig-prename descclassname\"><span class=\"pre\">structuralcodes.codes.mc2010.</span></span><span class=\"sig-name descname\"><span class=\"pre\">epsud</span></span><span class=\"sig-paren\">(</span><em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">epsuk</span></span><span class=\"p\"><span class=\"pre\">:</span></span><span class=\"w\"> </span><span class=\"n\"><span class=\"pre\">float</span></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">gamma_eps</span></span><span class=\"p\"><span class=\"pre\">:</span></span><span class=\"w\"> </span><span class=\"n\"><span class=\"pre\">float</span></span><span class=\"w\"> </span><span class=\"o\"><span class=\"pre\">=</span></span><span class=\"w\"> </span><span class=\"default_value\"><span class=\"pre\">0.9</span></span></em><span class=\"sig-paren\">)</span> <span class=\"sig-return\"><span class=\"sig-return-icon\">\u2192</span> <span class=\"sig-return-typehint\"><span class=\"pre\">float</span></span></span><a class=\"reference external\" href=\"https://github.com/fib-international/structuralcodes/tree/main/structuralcodes\\codes\\mc2010\\_reinforcement_material_properties.py#L48-L70\"><span class=\"viewcode-link\"><span class=\"pre\">[source]</span></span></a></dt><dd><p>Calculate the design value of the reinforcement ultimate strain.</p><p>fib Model Code 2010, Sec. 7.2.3.2.</p></dd>", "a[href=\"#material-properties-for-reinforcement-steel\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">Material properties for reinforcement steel<a class=\"headerlink\" href=\"#material-properties-for-reinforcement-steel\" title=\"Link to this heading\">\u00b6</a></h1>", "a[href=\"#structuralcodes.codes.mc2010.reinforcement_duct_props\"]": "<dt class=\"sig sig-object py\" id=\"structuralcodes.codes.mc2010.reinforcement_duct_props\">\n<span class=\"sig-prename descclassname\"><span class=\"pre\">structuralcodes.codes.mc2010.</span></span><span class=\"sig-name descname\"><span class=\"pre\">reinforcement_duct_props</span></span><span class=\"sig-paren\">(</span><em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">fyk</span></span><span class=\"p\"><span class=\"pre\">:</span></span><span class=\"w\"> </span><span class=\"n\"><span class=\"pre\">float</span></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">ductility_class</span></span><span class=\"p\"><span class=\"pre\">:</span></span><span class=\"w\"> </span><span class=\"n\"><span class=\"pre\">Literal</span><span class=\"p\"><span class=\"pre\">[</span></span><span class=\"s\"><span class=\"pre\">'A'</span></span><span class=\"p\"><span class=\"pre\">,</span></span><span class=\"w\"> </span><span class=\"s\"><span class=\"pre\">'B'</span></span><span class=\"p\"><span class=\"pre\">,</span></span><span class=\"w\"> </span><span class=\"s\"><span class=\"pre\">'C'</span></span><span class=\"p\"><span class=\"pre\">,</span></span><span class=\"w\"> </span><span class=\"s\"><span class=\"pre\">'D'</span></span><span class=\"p\"><span class=\"pre\">]</span></span></span></em><span class=\"sig-paren\">)</span> <span class=\"sig-return\"><span class=\"sig-return-icon\">\u2192</span> <span class=\"sig-return-typehint\"><span class=\"pre\">Dict</span><span class=\"p\"><span class=\"pre\">[</span></span><span class=\"pre\">str</span><span class=\"p\"><span class=\"pre\">,</span></span><span class=\"w\"> </span><span class=\"pre\">float</span><span class=\"p\"><span class=\"pre\">]</span></span></span></span><a class=\"reference external\" href=\"https://github.com/fib-international/structuralcodes/tree/main/structuralcodes\\codes\\mc2010\\_reinforcement_material_properties.py#L73-L105\"><span class=\"viewcode-link\"><span class=\"pre\">[source]</span></span></a></dt><dd><p>Return a dict with the minimum characteristic ductility properties for\nreinforcement ductility class.</p><p>fib Model Code 2010, Sec. 5.2.5.4.</p></dd>", "a[href=\"#structuralcodes.codes.mc2010.fyd\"]": "<dt class=\"sig sig-object py\" id=\"structuralcodes.codes.mc2010.fyd\">\n<span class=\"sig-prename descclassname\"><span class=\"pre\">structuralcodes.codes.mc2010.</span></span><span class=\"sig-name descname\"><span class=\"pre\">fyd</span></span><span class=\"sig-paren\">(</span><em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">fyk</span></span><span class=\"p\"><span class=\"pre\">:</span></span><span class=\"w\"> </span><span class=\"n\"><span class=\"pre\">float</span></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">gamma_s</span></span><span class=\"p\"><span class=\"pre\">:</span></span><span class=\"w\"> </span><span class=\"n\"><span class=\"pre\">float</span></span><span class=\"w\"> </span><span class=\"o\"><span class=\"pre\">=</span></span><span class=\"w\"> </span><span class=\"default_value\"><span class=\"pre\">1.15</span></span></em><span class=\"sig-paren\">)</span> <span class=\"sig-return\"><span class=\"sig-return-icon\">\u2192</span> <span class=\"sig-return-typehint\"><span class=\"pre\">float</span></span></span><a class=\"reference external\" href=\"https://github.com/fib-international/structuralcodes/tree/main/structuralcodes\\codes\\mc2010\\_reinforcement_material_properties.py#L25-L45\"><span class=\"viewcode-link\"><span class=\"pre\">[source]</span></span></a></dt><dd><p>Calculate the design value of the reinforcement yield strength.</p><p>fib Model Code 2010, Sec. 4.5.2.2.3.</p></dd>"}
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
