/*
 * TASSEL 5 nightly builds listing.
 *
 * Progressive enhancement: the page ships with static links to the rolling
 * `dev-latest` release plus a pointer at the GitHub releases page. When this
 * script can reach the GitHub Releases API it adds a table of every retained
 * nightly (the dated `dev-YYYYMMDD` prereleases). On any failure the static
 * content is left alone.
 */
(function () {
  "use strict";

  // Dated nightly tags only. The rolling `dev-latest` tag is a duplicate of the
  // newest of these and is covered by the static links above the table.
  var DATED_TAG = /^dev-(\d{4})(\d{2})(\d{2})$/;
  var ASSET = /^tassel-5-standalone-v(.+)-dev\.\d{8}\.(zip|tar\.gz)$/;
  var COMMIT_SHA = /^[0-9a-f]{40}$/i;

  var MONTHS = [
    "January", "February", "March", "April", "May", "June",
    "July", "August", "September", "October", "November", "December"
  ];

  function init() {
    var list = document.getElementById("nightly-list");
    if (!list) {
      return; // Not the nightly builds page.
    }

    var status = document.getElementById("nightly-status");
    var fallback = document.getElementById("nightly-fallback");

    if (!window.tasselReleases) {
      fail(list, status, new Error("releases.js did not load"));
      return;
    }

    window.tasselReleases
      .fetch()
      .then(function (releases) {
        var builds = collectBuilds(releases);
        if (builds.length === 0) {
          throw new Error("No nightly builds found");
        }
        list.setAttribute("data-state", "ready");
        list.innerHTML = renderTable(builds);
        if (fallback) {
          fallback.hidden = true;
        }
      })
      .catch(function (err) {
        fail(list, status, err);
      });
  }

  function fail(list, status, err) {
    if (window.console) {
      console.warn("TASSEL nightly list unavailable:", err);
    }
    list.setAttribute("data-state", "error");
    if (status) {
      status.textContent =
        "The live list of nightly builds is unavailable right now.";
    }
  }

  function collectBuilds(releases) {
    var builds = [];

    (releases || []).forEach(function (rel) {
      var match = DATED_TAG.exec(rel.tag_name || "");
      if (!rel.prerelease || !match) {
        return;
      }

      var build = {
        tag: rel.tag_name,
        date: match[1] + "-" + match[2] + "-" + match[3],
        label: MONTHS[parseInt(match[2], 10) - 1] + " " +
          parseInt(match[3], 10) + ", " + match[1],
        version: null,
        commit: COMMIT_SHA.test(rel.target_commitish || "")
          ? rel.target_commitish
          : null,
        assets: {}
      };

      (rel.assets || []).forEach(function (asset) {
        var m = ASSET.exec(asset.name || "");
        if (!m) {
          return;
        }
        build.version = build.version || m[1];
        build.assets[m[2] === "zip" ? "zip" : "targz"] = asset.browser_download_url;
      });

      if (build.assets.zip || build.assets.targz) {
        builds.push(build);
      }
    });

    builds.sort(function (a, b) {
      return a.date < b.date ? 1 : a.date > b.date ? -1 : 0;
    });

    return builds;
  }

  function renderTable(builds) {
    var repo = window.tasselReleases.repo;
    var rows = builds.map(function (build, index) {
      var links = [];
      if (build.assets.targz) {
        links.push(link(build.assets.targz, ".tar.gz"));
      }
      if (build.assets.zip) {
        links.push(link(build.assets.zip, ".zip"));
      }

      var commit = build.commit
        ? link(
            "https://github.com/" + repo + "/commit/" + build.commit,
            "<code>" + escapeHtml(build.commit.slice(0, 7)) + "</code>",
            true
          )
        : "&mdash;";

      return (
        "<tr>" +
        "<td>" + escapeHtml(build.label) +
        (index === 0 ? " <strong>(latest)</strong>" : "") + "</td>" +
        "<td>" + (build.version ? escapeHtml(build.version) : "&mdash;") + "</td>" +
        "<td>" + commit + "</td>" +
        "<td>" + links.join(" &middot; ") + "</td>" +
        "<td>" +
        link("https://github.com/" + repo + "/releases/tag/" + build.tag,
             "<code>" + escapeHtml(build.tag) + "</code>", true) +
        "</td>" +
        "</tr>"
      );
    });

    return (
      "<table>" +
      "<thead><tr>" +
      "<th>Built</th><th>Version</th><th>Commit</th>" +
      "<th>Standalone</th><th>Release</th>" +
      "</tr></thead>" +
      "<tbody>" + rows.join("") + "</tbody>" +
      "</table>"
    );
  }

  function link(url, label, isHtml) {
    return (
      '<a href="' + escapeAttr(url) + '">' +
      (isHtml ? label : escapeHtml(label)) +
      "</a>"
    );
  }

  function escapeHtml(s) {
    return String(s)
      .replace(/&/g, "&amp;")
      .replace(/</g, "&lt;")
      .replace(/>/g, "&gt;");
  }

  function escapeAttr(s) {
    return escapeHtml(s).replace(/"/g, "&quot;");
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", init);
  } else {
    init();
  }

  // Support Material for MkDocs instant navigation.
  if (window.document$ && typeof window.document$.subscribe === "function") {
    window.document$.subscribe(init);
  }
})();
