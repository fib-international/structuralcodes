selector_to_html = {"a[href=\"#structuralcodes.core.base.ConstitutiveLaw\"]": "<dt class=\"sig sig-object py\" id=\"structuralcodes.core.base.ConstitutiveLaw\">\n<em class=\"property\"><span class=\"pre\">class</span><span class=\"w\"> </span></em><span class=\"sig-prename descclassname\"><span class=\"pre\">structuralcodes.core.base.</span></span><span class=\"sig-name descname\"><span class=\"pre\">ConstitutiveLaw</span></span><span class=\"sig-paren\">(</span><em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">name</span></span><span class=\"p\"><span class=\"pre\">:</span></span><span class=\"w\"> </span><span class=\"n\"><span class=\"pre\">str</span><span class=\"w\"> </span><span class=\"p\"><span class=\"pre\">|</span></span><span class=\"w\"> </span><span class=\"pre\">None</span></span><span class=\"w\"> </span><span class=\"o\"><span class=\"pre\">=</span></span><span class=\"w\"> </span><span class=\"default_value\"><span class=\"pre\">None</span></span></em><span class=\"sig-paren\">)</span><a class=\"reference external\" href=\"https://github.com/fib-international/structuralcodes/tree/main/structuralcodes\\core\\base.py#L68-L192\"><span class=\"viewcode-link\"><span class=\"pre\">[source]</span></span></a></dt><dd><p>Abstract base class for constitutive laws.</p></dd>", "a[href=\"#structuralcodes.core.base.Material.update_attributes\"]": "<dt class=\"sig sig-object py\" id=\"structuralcodes.core.base.Material.update_attributes\">\n<span class=\"sig-name descname\"><span class=\"pre\">update_attributes</span></span><span class=\"sig-paren\">(</span><em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">updated_attributes</span></span><span class=\"p\"><span class=\"pre\">:</span></span><span class=\"w\"> </span><span class=\"n\"><span class=\"pre\">Dict</span></span></em><span class=\"sig-paren\">)</span> <span class=\"sig-return\"><span class=\"sig-return-icon\">\u2192</span> <span class=\"sig-return-typehint\"><span class=\"pre\">None</span></span></span><a class=\"reference external\" href=\"https://github.com/fib-international/structuralcodes/tree/main/structuralcodes\\core\\base.py#L32-L50\"><span class=\"viewcode-link\"><span class=\"pre\">[source]</span></span></a></dt><dd><p>Function for updating the attributes specified in the input\ndictionary.</p></dd>", "a[href=\"#structuralcodes.materials.reinforcement.Reinforcement.epsyd\"]": "<dt class=\"sig sig-object py\" id=\"structuralcodes.materials.reinforcement.Reinforcement.epsyd\">\n<em class=\"property\"><span class=\"pre\">property</span><span class=\"w\"> </span></em><span class=\"sig-name descname\"><span class=\"pre\">epsyd</span></span><em class=\"property\"><span class=\"p\"><span class=\"pre\">:</span></span><span class=\"w\"> </span><span class=\"pre\">float</span></em></dt><dd><p>Return the design yield strain.</p></dd>", "a[href=\"#structuralcodes.core.base.Material\"]": "<dt class=\"sig sig-object py\" id=\"structuralcodes.core.base.Material\">\n<em class=\"property\"><span class=\"pre\">class</span><span class=\"w\"> </span></em><span class=\"sig-prename descclassname\"><span class=\"pre\">structuralcodes.core.base.</span></span><span class=\"sig-name descname\"><span class=\"pre\">Material</span></span><span class=\"sig-paren\">(</span><em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">density</span></span><span class=\"p\"><span class=\"pre\">:</span></span><span class=\"w\"> </span><span class=\"n\"><span class=\"pre\">float</span></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">name</span></span><span class=\"p\"><span class=\"pre\">:</span></span><span class=\"w\"> </span><span class=\"n\"><span class=\"pre\">str</span><span class=\"w\"> </span><span class=\"p\"><span class=\"pre\">|</span></span><span class=\"w\"> </span><span class=\"pre\">None</span></span><span class=\"w\"> </span><span class=\"o\"><span class=\"pre\">=</span></span><span class=\"w\"> </span><span class=\"default_value\"><span class=\"pre\">None</span></span></em><span class=\"sig-paren\">)</span><a class=\"reference external\" href=\"https://github.com/fib-international/structuralcodes/tree/main/structuralcodes\\core\\base.py#L15-L65\"><span class=\"viewcode-link\"><span class=\"pre\">[source]</span></span></a></dt><dd><p>Abstract base class for materials.</p></dd>", "a[href=\"#structuralcodes.materials.concrete.Concrete\"]": "<dt class=\"sig sig-object py\" id=\"structuralcodes.materials.concrete.Concrete\">\n<em class=\"property\"><span class=\"pre\">class</span><span class=\"w\"> </span></em><span class=\"sig-prename descclassname\"><span class=\"pre\">structuralcodes.materials.concrete.</span></span><span class=\"sig-name descname\"><span class=\"pre\">Concrete</span></span><span class=\"sig-paren\">(</span><em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">fck</span></span><span class=\"p\"><span class=\"pre\">:</span></span><span class=\"w\"> </span><span class=\"n\"><span class=\"pre\">float</span></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">name</span></span><span class=\"p\"><span class=\"pre\">:</span></span><span class=\"w\"> </span><span class=\"n\"><span class=\"pre\">str</span><span class=\"w\"> </span><span class=\"p\"><span class=\"pre\">|</span></span><span class=\"w\"> </span><span class=\"pre\">None</span></span><span class=\"w\"> </span><span class=\"o\"><span class=\"pre\">=</span></span><span class=\"w\"> </span><span class=\"default_value\"><span class=\"pre\">None</span></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">density</span></span><span class=\"p\"><span class=\"pre\">:</span></span><span class=\"w\"> </span><span class=\"n\"><span class=\"pre\">float</span></span><span class=\"w\"> </span><span class=\"o\"><span class=\"pre\">=</span></span><span class=\"w\"> </span><span class=\"default_value\"><span class=\"pre\">2400</span></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">gamma_c</span></span><span class=\"p\"><span class=\"pre\">:</span></span><span class=\"w\"> </span><span class=\"n\"><span class=\"pre\">float</span><span class=\"w\"> </span><span class=\"p\"><span class=\"pre\">|</span></span><span class=\"w\"> </span><span class=\"pre\">None</span></span><span class=\"w\"> </span><span class=\"o\"><span class=\"pre\">=</span></span><span class=\"w\"> </span><span class=\"default_value\"><span class=\"pre\">None</span></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">existing</span></span><span class=\"p\"><span class=\"pre\">:</span></span><span class=\"w\"> </span><span class=\"n\"><span class=\"pre\">bool</span><span class=\"w\"> </span><span class=\"p\"><span class=\"pre\">|</span></span><span class=\"w\"> </span><span class=\"pre\">None</span></span><span class=\"w\"> </span><span class=\"o\"><span class=\"pre\">=</span></span><span class=\"w\"> </span><span class=\"default_value\"><span class=\"pre\">False</span></span></em><span class=\"sig-paren\">)</span><a class=\"reference external\" href=\"https://github.com/fib-international/structuralcodes/tree/main/structuralcodes\\materials\\concrete\\_concrete.py#L10-L114\"><span class=\"viewcode-link\"><span class=\"pre\">[source]</span></span></a></dt><dd><p>The abstract concrete material.</p></dd>", "a[href=\"#structuralcodes.core.base.ConstitutiveLaw.get_secant\"]": "<dt class=\"sig sig-object py\" id=\"structuralcodes.core.base.ConstitutiveLaw.get_secant\">\n<span class=\"sig-name descname\"><span class=\"pre\">get_secant</span></span><span class=\"sig-paren\">(</span><em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">eps</span></span><span class=\"p\"><span class=\"pre\">:</span></span><span class=\"w\"> </span><span class=\"n\"><span class=\"pre\">float</span></span></em><span class=\"sig-paren\">)</span> <span class=\"sig-return\"><span class=\"sig-return-icon\">\u2192</span> <span class=\"sig-return-typehint\"><span class=\"pre\">float</span></span></span><a class=\"reference external\" href=\"https://github.com/fib-international/structuralcodes/tree/main/structuralcodes\\core\\base.py#L185-L192\"><span class=\"viewcode-link\"><span class=\"pre\">[source]</span></span></a></dt><dd><p>Method to return the\nsecant at a given strain level.</p></dd>", "a[href=\"#structuralcodes.materials.concrete.Concrete.constitutive_law\"]": "<dt class=\"sig sig-object py\" id=\"structuralcodes.materials.concrete.Concrete.constitutive_law\">\n<em class=\"property\"><span class=\"pre\">property</span><span class=\"w\"> </span></em><span class=\"sig-name descname\"><span class=\"pre\">constitutive_law</span></span><em class=\"property\"><span class=\"p\"><span class=\"pre\">:</span></span><span class=\"w\"> </span><a class=\"reference internal\" href=\"#structuralcodes.core.base.ConstitutiveLaw\" title=\"structuralcodes.core.base.ConstitutiveLaw\"><span class=\"pre\">ConstitutiveLaw</span></a></em></dt><dd><p>Returns the constitutive law object.</p></dd>", "a[href=\"#structuralcodes.materials.reinforcement.Reinforcement\"]": "<dt class=\"sig sig-object py\" id=\"structuralcodes.materials.reinforcement.Reinforcement\">\n<em class=\"property\"><span class=\"pre\">class</span><span class=\"w\"> </span></em><span class=\"sig-prename descclassname\"><span class=\"pre\">structuralcodes.materials.reinforcement.</span></span><span class=\"sig-name descname\"><span class=\"pre\">Reinforcement</span></span><span class=\"sig-paren\">(</span><em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">fyk</span></span><span class=\"p\"><span class=\"pre\">:</span></span><span class=\"w\"> </span><span class=\"n\"><span class=\"pre\">float</span></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">Es</span></span><span class=\"p\"><span class=\"pre\">:</span></span><span class=\"w\"> </span><span class=\"n\"><span class=\"pre\">float</span></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">density</span></span><span class=\"p\"><span class=\"pre\">:</span></span><span class=\"w\"> </span><span class=\"n\"><span class=\"pre\">float</span></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">ftk</span></span><span class=\"p\"><span class=\"pre\">:</span></span><span class=\"w\"> </span><span class=\"n\"><span class=\"pre\">float</span></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">epsuk</span></span><span class=\"p\"><span class=\"pre\">:</span></span><span class=\"w\"> </span><span class=\"n\"><span class=\"pre\">float</span></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">gamma_s</span></span><span class=\"p\"><span class=\"pre\">:</span></span><span class=\"w\"> </span><span class=\"n\"><span class=\"pre\">float</span><span class=\"w\"> </span><span class=\"p\"><span class=\"pre\">|</span></span><span class=\"w\"> </span><span class=\"pre\">None</span></span><span class=\"w\"> </span><span class=\"o\"><span class=\"pre\">=</span></span><span class=\"w\"> </span><span class=\"default_value\"><span class=\"pre\">None</span></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">name</span></span><span class=\"p\"><span class=\"pre\">:</span></span><span class=\"w\"> </span><span class=\"n\"><span class=\"pre\">str</span><span class=\"w\"> </span><span class=\"p\"><span class=\"pre\">|</span></span><span class=\"w\"> </span><span class=\"pre\">None</span></span><span class=\"w\"> </span><span class=\"o\"><span class=\"pre\">=</span></span><span class=\"w\"> </span><span class=\"default_value\"><span class=\"pre\">None</span></span></em><span class=\"sig-paren\">)</span><a class=\"reference external\" href=\"https://github.com/fib-international/structuralcodes/tree/main/structuralcodes\\materials\\reinforcement\\_reinforcement.py#L13-L172\"><span class=\"viewcode-link\"><span class=\"pre\">[source]</span></span></a></dt><dd><p>The abstract reinforcement material.</p></dd>", "a[href=\"#base-concrete-class\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Base concrete class<a class=\"headerlink\" href=\"#base-concrete-class\" title=\"Link to this heading\">\u00b6</a></h2>", "a[href=\"#base-material-class\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Base material class<a class=\"headerlink\" href=\"#base-material-class\" title=\"Link to this heading\">\u00b6</a></h2>", "a[href=\"#structuralcodes.materials.reinforcement.Reinforcement.epsyk\"]": "<dt class=\"sig sig-object py\" id=\"structuralcodes.materials.reinforcement.Reinforcement.epsyk\">\n<em class=\"property\"><span class=\"pre\">property</span><span class=\"w\"> </span></em><span class=\"sig-name descname\"><span class=\"pre\">epsyk</span></span><em class=\"property\"><span class=\"p\"><span class=\"pre\">:</span></span><span class=\"w\"> </span><span class=\"pre\">float</span></em></dt><dd><p>Returns characteristic yield strain epsyk.</p></dd>", "a[href=\"#base-constitutive-law-class\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Base constitutive law class<a class=\"headerlink\" href=\"#base-constitutive-law-class\" title=\"Link to this heading\">\u00b6</a></h2>", "a[href=\"#structuralcodes.core.base.ConstitutiveLaw.name\"]": "<dt class=\"sig sig-object py\" id=\"structuralcodes.core.base.ConstitutiveLaw.name\">\n<em class=\"property\"><span class=\"pre\">property</span><span class=\"w\"> </span></em><span class=\"sig-name descname\"><span class=\"pre\">name</span></span></dt><dd><p>Returns the name of the constitutive law.</p></dd>", "a[href=\"#structuralcodes.materials.reinforcement.Reinforcement.constitutive_law\"]": "<dt class=\"sig sig-object py\" id=\"structuralcodes.materials.reinforcement.Reinforcement.constitutive_law\">\n<em class=\"property\"><span class=\"pre\">property</span><span class=\"w\"> </span></em><span class=\"sig-name descname\"><span class=\"pre\">constitutive_law</span></span><em class=\"property\"><span class=\"p\"><span class=\"pre\">:</span></span><span class=\"w\"> </span><a class=\"reference internal\" href=\"#structuralcodes.core.base.ConstitutiveLaw\" title=\"structuralcodes.core.base.ConstitutiveLaw\"><span class=\"pre\">ConstitutiveLaw</span></a></em></dt><dd><p>Returns the constitutive law object.</p></dd>", "a[href=\"#base-material-classes\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">Base material classes<a class=\"headerlink\" href=\"#base-material-classes\" title=\"Link to this heading\">\u00b6</a></h1><h2>Base material class<a class=\"headerlink\" href=\"#base-material-class\" title=\"Link to this heading\">\u00b6</a></h2>", "a[href=\"#structuralcodes.materials.reinforcement.Reinforcement.ftk\"]": "<dt class=\"sig sig-object py\" id=\"structuralcodes.materials.reinforcement.Reinforcement.ftk\">\n<em class=\"property\"><span class=\"pre\">property</span><span class=\"w\"> </span></em><span class=\"sig-name descname\"><span class=\"pre\">ftk</span></span><em class=\"property\"><span class=\"p\"><span class=\"pre\">:</span></span><span class=\"w\"> </span><span class=\"pre\">float</span></em></dt><dd><p>Returns ftk in MPa.</p></dd>", "a[href=\"#structuralcodes.materials.reinforcement.Reinforcement.Es\"]": "<dt class=\"sig sig-object py\" id=\"structuralcodes.materials.reinforcement.Reinforcement.Es\">\n<em class=\"property\"><span class=\"pre\">property</span><span class=\"w\"> </span></em><span class=\"sig-name descname\"><span class=\"pre\">Es</span></span><em class=\"property\"><span class=\"p\"><span class=\"pre\">:</span></span><span class=\"w\"> </span><span class=\"pre\">float</span></em></dt><dd><p>Returns Es in MPa.</p></dd>", "a[href=\"#structuralcodes.core.base.Material.density\"]": "<dt class=\"sig sig-object py\" id=\"structuralcodes.core.base.Material.density\">\n<em class=\"property\"><span class=\"pre\">property</span><span class=\"w\"> </span></em><span class=\"sig-name descname\"><span class=\"pre\">density</span></span></dt><dd><p>Returns the density of the material in kg/m3.</p></dd>", "a[href=\"#structuralcodes.materials.reinforcement.Reinforcement.fyk\"]": "<dt class=\"sig sig-object py\" id=\"structuralcodes.materials.reinforcement.Reinforcement.fyk\">\n<em class=\"property\"><span class=\"pre\">property</span><span class=\"w\"> </span></em><span class=\"sig-name descname\"><span class=\"pre\">fyk</span></span><em class=\"property\"><span class=\"p\"><span class=\"pre\">:</span></span><span class=\"w\"> </span><span class=\"pre\">float</span></em></dt><dd><p>Returns fyk in MPa.</p></dd>", "a[href=\"#structuralcodes.core.base.ConstitutiveLaw.preprocess_strains_with_limits\"]": "<dt class=\"sig sig-object py\" id=\"structuralcodes.core.base.ConstitutiveLaw.preprocess_strains_with_limits\">\n<span class=\"sig-name descname\"><span class=\"pre\">preprocess_strains_with_limits</span></span><span class=\"sig-paren\">(</span><em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">eps</span></span><span class=\"p\"><span class=\"pre\">:</span></span><span class=\"w\"> </span><span class=\"n\"><span class=\"pre\">float</span><span class=\"w\"> </span><span class=\"p\"><span class=\"pre\">|</span></span><span class=\"w\"> </span><span class=\"pre\">ArrayLike</span></span></em><span class=\"sig-paren\">)</span> <span class=\"sig-return\"><span class=\"sig-return-icon\">\u2192</span> <span class=\"sig-return-typehint\"><span class=\"pre\">float</span><span class=\"w\"> </span><span class=\"p\"><span class=\"pre\">|</span></span><span class=\"w\"> </span><span class=\"pre\">ArrayLike</span></span></span><a class=\"reference external\" href=\"https://github.com/fib-international/structuralcodes/tree/main/structuralcodes\\core\\base.py#L111-L131\"><span class=\"viewcode-link\"><span class=\"pre\">[source]</span></span></a></dt><dd><p>Preprocess strain arrays setting those strains sufficiently\nnear to ultimate strain limits to exactly ultimate strain limit.</p></dd>", "a[href=\"#base-reinforcement-steel-class\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Base reinforcement steel class<a class=\"headerlink\" href=\"#base-reinforcement-steel-class\" title=\"Link to this heading\">\u00b6</a></h2>", "a[href=\"#structuralcodes.materials.concrete.Concrete.fck\"]": "<dt class=\"sig sig-object py\" id=\"structuralcodes.materials.concrete.Concrete.fck\">\n<em class=\"property\"><span class=\"pre\">property</span><span class=\"w\"> </span></em><span class=\"sig-name descname\"><span class=\"pre\">fck</span></span><em class=\"property\"><span class=\"p\"><span class=\"pre\">:</span></span><span class=\"w\"> </span><span class=\"pre\">float</span></em></dt><dd><p>Returns fck in MPa.</p></dd>", "a[href=\"#structuralcodes.core.base.Material.name\"]": "<dt class=\"sig sig-object py\" id=\"structuralcodes.core.base.Material.name\">\n<em class=\"property\"><span class=\"pre\">property</span><span class=\"w\"> </span></em><span class=\"sig-name descname\"><span class=\"pre\">name</span></span></dt><dd><p>Returns the name of the material.</p></dd>", "a[href=\"#structuralcodes.materials.concrete.Concrete.gamma_c\"]": "<dt class=\"sig sig-object py\" id=\"structuralcodes.materials.concrete.Concrete.gamma_c\">\n<em class=\"property\"><span class=\"pre\">abstract</span><span class=\"w\"> </span><span class=\"pre\">property</span><span class=\"w\"> </span></em><span class=\"sig-name descname\"><span class=\"pre\">gamma_c</span></span><em class=\"property\"><span class=\"p\"><span class=\"pre\">:</span></span><span class=\"w\"> </span><span class=\"pre\">float</span></em></dt><dd><p>Each concrete should implement its own getter for the partial factor\nin order to interact with the globally set national annex.</p></dd>", "a[href=\"#structuralcodes.materials.reinforcement.Reinforcement.epsuk\"]": "<dt class=\"sig sig-object py\" id=\"structuralcodes.materials.reinforcement.Reinforcement.epsuk\">\n<em class=\"property\"><span class=\"pre\">property</span><span class=\"w\"> </span></em><span class=\"sig-name descname\"><span class=\"pre\">epsuk</span></span><em class=\"property\"><span class=\"p\"><span class=\"pre\">:</span></span><span class=\"w\"> </span><span class=\"pre\">float</span></em></dt><dd><p>Returns epsuk.</p></dd>", "a[href=\"#structuralcodes.core.base.Material.constitutive_law\"]": "<dt class=\"sig sig-object py\" id=\"structuralcodes.core.base.Material.constitutive_law\">\n<em class=\"property\"><span class=\"pre\">property</span><span class=\"w\"> </span></em><span class=\"sig-name descname\"><span class=\"pre\">constitutive_law</span></span></dt><dd><p>Returns the ConstitutiveLaw of the object.</p></dd>"}
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