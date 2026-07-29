/*
 * Style changelog version tags in section headings as monospace.
 */
(function () {
  "use strict";

  function escapeHtml(s) {
    return String(s)
      .replace(/&/g, "&amp;")
      .replace(/</g, "&lt;")
      .replace(/>/g, "&gt;");
  }

  function cleanHeadingText(text) {
    return String(text)
      .replace(/\u00B6/g, "")
      .replace(/¶/g, "")
      .trim();
  }

  function styleVersionTags() {
    if (!/\/changelog\/?$/.test(location.pathname)) {
      return;
    }

    document.body.classList.add("page-changelog");

    document.querySelectorAll(".md-typeset h2").forEach(function (h2) {
      h2.querySelectorAll(".headerlink").forEach(function (link) {
        link.remove();
      });

      var text = cleanHeadingText(h2.textContent || "");
      var match = text.match(/^\(V([^)]+)\)\s+(.*)$/);
      if (!match) {
        return;
      }

      var tagIcon =
        '<span class="material-symbols-outlined changelog-version__icon" aria-hidden="true">sell</span>';

      h2.innerHTML =
        '<span class="changelog-version">' +
        tagIcon +
        '<span class="changelog-version__text">v' +
        escapeHtml(match[1]) +
        "</span></span> " +
        '<span class="changelog-date">' +
        escapeHtml(cleanHeadingText(match[2])) +
        "</span>";
      h2.dataset.versionStyled = "true";
    });
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", styleVersionTags);
  } else {
    styleVersionTags();
  }

  if (window.document$ && typeof window.document$.subscribe === "function") {
    window.document$.subscribe(styleVersionTags);
  }
})();
