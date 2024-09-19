selector_to_html = {"a[href=\"#set-up-your-local-development-environment\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Set up your local development environment<a class=\"headerlink\" href=\"#set-up-your-local-development-environment\" title=\"Link to this heading\">\u00b6</a></h2><p>To keep track of the library dependencies, we recommend to work in a virtual environment. To create one on Windows, type <code class=\"docutils literal notranslate\"><span class=\"pre\">py</span> <span class=\"pre\">-m</span> <span class=\"pre\">venv</span> <span class=\"pre\">venv</span></code>, and activate it by typing <code class=\"docutils literal notranslate\"><span class=\"pre\">venv\\scripts\\activate</span></code>.</p><p>We use <code class=\"docutils literal notranslate\"><span class=\"pre\">ruff</span></code> for linting and code formatting, and have set up our test suite using <code class=\"docutils literal notranslate\"><span class=\"pre\">pytest</span></code>. Make sure you have these installed (<code class=\"docutils literal notranslate\"><span class=\"pre\">make</span> <span class=\"pre\">deps</span></code>). We have set up VS Code to use testing. Open VS Code inside the repo (<code class=\"docutils literal notranslate\"><span class=\"pre\">code</span> <span class=\"pre\">.</span></code>) and view the tests in the test pane on the left side.</p>", "a[href=\"#contribution-guidelines\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Contribution guidelines<a class=\"headerlink\" href=\"#contribution-guidelines\" title=\"Link to this heading\">\u00b6</a></h2><p>Here is how to get started.</p><p><em>First</em></p>", "a[href=\"#contributing\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">Contributing<a class=\"headerlink\" href=\"#contributing\" title=\"Link to this heading\">\u00b6</a></h1><p><a class=\"reference external\" href=\"https://github.com/fib-international/structuralcodes/actions/workflows/build.yaml\"><img alt=\"Build\" src=\"https://github.com/fib-international/structuralcodes/actions/workflows/build.yaml/badge.svg\"/></a>\n<a class=\"reference external\" href=\"https://github.com/charliermarsh/ruff\"><img alt=\"Code style: Ruff\" src=\"https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/charliermarsh/ruff/main/assets/badge/format.json\"/></a>\n<a class=\"reference external\" href=\"https://github.com/charliermarsh/ruff\"><img alt=\"Ruff\" src=\"https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/charliermarsh/ruff/main/assets/badge/v2.json\"/></a></p>", "a[href=\"#local-dev-env\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Set up your local development environment<a class=\"headerlink\" href=\"#set-up-your-local-development-environment\" title=\"Link to this heading\">\u00b6</a></h2><p>To keep track of the library dependencies, we recommend to work in a virtual environment. To create one on Windows, type <code class=\"docutils literal notranslate\"><span class=\"pre\">py</span> <span class=\"pre\">-m</span> <span class=\"pre\">venv</span> <span class=\"pre\">venv</span></code>, and activate it by typing <code class=\"docutils literal notranslate\"><span class=\"pre\">venv\\scripts\\activate</span></code>.</p><p>We use <code class=\"docutils literal notranslate\"><span class=\"pre\">ruff</span></code> for linting and code formatting, and have set up our test suite using <code class=\"docutils literal notranslate\"><span class=\"pre\">pytest</span></code>. Make sure you have these installed (<code class=\"docutils literal notranslate\"><span class=\"pre\">make</span> <span class=\"pre\">deps</span></code>). We have set up VS Code to use testing. Open VS Code inside the repo (<code class=\"docutils literal notranslate\"><span class=\"pre\">code</span> <span class=\"pre\">.</span></code>) and view the tests in the test pane on the left side.</p>", "a[href=\"#general\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">General<a class=\"headerlink\" href=\"#general\" title=\"Link to this heading\">\u00b6</a></h2><p>For the moment, the core team has completed an initial structure for the two lowest levels, i.e. how to structure the models from the structural design codes, and how to create classes that represent materials and sections. Contributions that seek to add content to the lowest level are highly appreciated.</p><p>We have defined a set of <a class=\"reference external\" href=\"https://github.com/fib-international/structuralcodes/milestones\">milestones</a> that indicate our priorities. Issues are linked to each milesstone, and some of the issues include a to-do list and have linked pull requests. Contributions related to the milestones are highly appreciated. For the moment, our first priority is to add code related to <em>fib</em> Model Code 2020, but other contributions are of course also welcome!</p>"}
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
