// sphinx_design renders the `dropdown` directive as a <details> element,
// which browsers do not automatically expand when navigated to via a URL
// fragment (e.g., from the parameter-set class index above). This opens the
// target (or its nearest <details> ancestor) on load and on same-page hash
// navigation, then re-scrolls it into view since opening the dropdown shifts
// the page layout.
//
// This file is only loaded by doc/pypeit_par.rst (see build_par_rst.py); it
// is intentionally not registered site-wide in conf.py's html_js_files, so
// this behavior is confined to the parameter-set tables page.
document.addEventListener('DOMContentLoaded', function () {
    function openTargetDropdown() {
        if (!window.location.hash) {
            return;
        }
        var target;
        try {
            target = document.querySelector(window.location.hash);
        } catch (e) {
            return;
        }
        if (!target) {
            return;
        }
        var details = target.closest('details');
        if (details && !details.open) {
            details.open = true;
            target.scrollIntoView();
        }
    }
    openTargetDropdown();
    window.addEventListener('hashchange', openTargetDropdown);
});
