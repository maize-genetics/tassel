/*
 * TASSEL 5 download page enhancement.
 *
 * Progressive enhancement: the page ships with a static fallback list of direct
 * download links. When this script can reach the GitHub Releases API it replaces
 * the static list with a live, auto-detecting picker (OS detection + full
 * version history). On any failure it silently leaves the static fallback in
 * place.
 */
(function () {
  "use strict";

  var REPO = "maize-genetics/tassel";
  var API = "https://api.github.com/repos/" + REPO + "/releases?per_page=100";
  var RELEASES_URL = "https://github.com/" + REPO + "/releases";
  var CACHE_KEY = "tassel-releases-v1";
  var CACHE_TTL_MS = 30 * 60 * 1000; // 30 minutes

  // Installer OS options (native jDeploy installers are attached to each version release).
  var INSTALLER_OS = [
    { value: "mac-arm64", label: "macOS (Apple Silicon)", match: "mac-arm64" },
    { value: "mac-x64", label: "macOS (Intel)", match: "mac-x64" },
    { value: "linux-x64", label: "Linux / Unix", match: "linux-x64" },
    { value: "win-x64", label: "Windows (64-bit)", match: "win-x64" }
  ];

  function $(id) {
    return document.getElementById(id);
  }

  function init() {
    var distSel = $("dl-dist");
    if (!distSel) {
      return; // Not the download page.
    }

    var els = {
      dist: distSel,
      os: $("dl-os"),
      osField: $("dl-os-field"),
      version: $("dl-version"),
      result: $("dl-result"),
      fallback: $("dl-fallback"),
      recWrap: $("dl-recommended"),
      recTitle: $("dl-rec-title"),
      recMeta: $("dl-rec-meta"),
      recBtn: $("dl-rec-btn"),
      cliWrap: $("dl-standalone"),
      cliTitle: $("dl-cli-title"),
      cliMeta: $("dl-cli-meta"),
      cliBtn: $("dl-cli-btn")
    };

    fetchReleases()
      .then(function (releases) {
        var data = buildCatalog(releases);
        if (!data.installers && data.standaloneVersions.length === 0) {
          throw new Error("No downloadable assets found");
        }
        observeFileTicker(els.recMeta);
        observeFileTicker(els.cliMeta);
        wirePicker(els, data);
        setRecommended(els, data);
        setStandaloneRecommended(els, data);
        if (els.fallback) {
          els.fallback.hidden = true;
        }
      })
      .catch(function (err) {
        // Keep the static fallback visible.
        if (window.console) {
          console.warn("TASSEL download picker unavailable:", err);
        }
        if (els.recWrap) {
          els.recWrap.setAttribute("data-state", "error");
          els.recTitle.textContent = "Desktop application";
          els.recMeta.textContent =
            "Automatic detection is unavailable right now.";
          els.recBtn.setAttribute("href", "#dl-fallback");
          els.recBtn.textContent = "See direct downloads";
        }
        if (els.cliWrap) {
          els.cliWrap.setAttribute("data-state", "error");
          els.cliTitle.textContent = "Standalone CLI";
          els.cliMeta.textContent =
            "Automatic detection is unavailable right now.";
          els.cliBtn.setAttribute("href", "#dl-fallback");
          els.cliBtn.textContent = "See direct downloads";
        }
        var picker = document.querySelector(".dl-picker");
        if (picker) {
          picker.hidden = true;
        }
      });
  }

  function fetchReleases() {
    try {
      var cached = sessionStorage.getItem(CACHE_KEY);
      if (cached) {
        var parsed = JSON.parse(cached);
        if (parsed && Date.now() - parsed.t < CACHE_TTL_MS) {
          return Promise.resolve(parsed.d);
        }
      }
    } catch (e) {
      /* ignore cache errors */
    }

    return fetch(API, { headers: { Accept: "application/vnd.github+json" } })
      .then(function (resp) {
        if (!resp.ok) {
          throw new Error("GitHub API " + resp.status);
        }
        return resp.json();
      })
      .then(function (releases) {
        try {
          sessionStorage.setItem(
            CACHE_KEY,
            JSON.stringify({ t: Date.now(), d: releases })
          );
        } catch (e) {
          /* ignore */
        }
        return releases;
      });
  }

  function buildCatalog(releases) {
    var installers = null; // { "mac-arm64": {name, url}, ... }
    var installerTag = null; // Release tag the chosen installers came from.
    var standalone = {}; // version -> { targz, zip }

    (releases || []).forEach(function (rel) {
      var assets = rel.assets || [];
      assets.forEach(function (asset) {
        var name = asset.name || "";
        var url = asset.browser_download_url;

        // Native installers, attached to each version release.
        if (name.indexOf("TASSEL.5.Installer-") === 0) {
          INSTALLER_OS.forEach(function (os) {
            if (name.indexOf(os.match) !== -1) {
              installers = installers || {};
              // The API lists releases newest first, so the first match wins.
              if (!installers[os.value]) {
                installers[os.value] = { name: name, url: url };
                installerTag = installerTag || rel.tag_name;
              }
            }
          });
        }

        // Standalone distributions (per version tag).
        var m = name.match(/tassel-5-standalone-v([0-9][0-9A-Za-z.\-]*)\.(tar\.gz|zip)$/);
        if (m) {
          var ver = m[1];
          var kind = m[2] === "zip" ? "zip" : "targz";
          standalone[ver] = standalone[ver] || {};
          standalone[ver][kind] = { name: name, url: url };
        }
      });
    });

    var standaloneVersions = Object.keys(standalone).sort(compareVersionsDesc);

    return {
      installers: installers,
      installerTag: installerTag,
      standalone: standalone,
      standaloneVersions: standaloneVersions
    };
  }

  function compareVersionsDesc(a, b) {
    var pa = a.split(/[.\-]/).map(numOr);
    var pb = b.split(/[.\-]/).map(numOr);
    for (var i = 0; i < Math.max(pa.length, pb.length); i++) {
      var da = pa[i] || 0;
      var db = pb[i] || 0;
      if (da !== db) {
        return db - da;
      }
    }
    return 0;
  }

  function numOr(x) {
    var n = parseInt(x, 10);
    return isNaN(n) ? 0 : n;
  }

  function wirePicker(els, data) {
    function onDistChange() {
      var dist = els.dist.value;
      if (dist === "installer") {
        els.osField.hidden = false;
        fillSelect(
          els.os,
          INSTALLER_OS.filter(function (o) {
            return data.installers && data.installers[o.value];
          }).map(function (o) {
            return { value: o.value, label: o.label };
          })
        );
        fillSelect(els.version, [
          {
            value: data.installerTag || "latest",
            label: data.installerTag
              ? "Latest (" + data.installerTag + ")"
              : "Latest"
          }
        ]);
      } else {
        els.osField.hidden = true;
        fillSelect(els.version,
          data.standaloneVersions.map(function (v) {
            return { value: v, label: "v" + v };
          })
        );
      }
      render();
    }

    function render() {
      var dist = els.dist.value;
      if (dist === "installer") {
        renderInstaller(els, data, els.os.value);
      } else {
        renderStandalone(els, data, els.version.value);
      }
    }

    els.dist.addEventListener("change", onDistChange);
    els.os.addEventListener("change", render);
    els.version.addEventListener("change", render);

    onDistChange();
  }

  function renderInstaller(els, data, osValue) {
    var asset = data.installers && data.installers[osValue];
    var osLabel = labelFor(INSTALLER_OS, osValue);
    if (!asset) {
      els.result.innerHTML =
        '<p class="dl-result__empty">No native installer is available for ' +
        escapeHtml(osLabel) +
        ". Try the standalone distribution.</p>";
      return;
    }
    els.result.innerHTML = downloadCard({
      title: "TASSEL 5 for " + osLabel,
      filename: asset.name,
      url: asset.url,
      note:
        "Native installer &middot; " +
        (data.installerTag ? escapeHtml(data.installerTag) : "latest build")
    });
  }

  function renderStandalone(els, data, version) {
    var entry = data.standalone[version];
    if (!entry) {
      els.result.innerHTML =
        '<p class="dl-result__empty">No standalone build found for this version.</p>';
      return;
    }
    var buttons = "";
    if (entry.targz) {
      buttons += buttonLink(entry.targz.url, "Download .tar.gz");
    }
    if (entry.zip) {
      buttons += buttonLink(entry.zip.url, "Download .zip");
    }
    var primary = entry.targz || entry.zip;
    els.result.innerHTML =
      '<div class="dl-card">' +
      '<div class="dl-card__body">' +
      '<h3 class="dl-card__title">TASSEL 5 Standalone v' +
      escapeHtml(version) +
      "</h3>" +
      '<p class="dl-card__meta">Command-line distribution &middot; requires Java 21+</p>' +
      '<p class="dl-card__file">' + escapeHtml(primary.name) + "</p>" +
      "</div>" +
      '<div class="dl-card__actions">' + buttons + "</div>" +
      "</div>";
  }

  function downloadCard(opts) {
    return (
      '<div class="dl-card">' +
      '<div class="dl-card__body">' +
      '<h3 class="dl-card__title">' + escapeHtml(opts.title) + "</h3>" +
      '<p class="dl-card__meta">' + opts.note + "</p>" +
      '<p class="dl-card__file">' + escapeHtml(opts.filename) + "</p>" +
      "</div>" +
      '<div class="dl-card__actions">' +
      '<a class="md-button md-button--primary" href="' +
      escapeAttr(opts.url) +
      '">Download</a>' +
      "</div>" +
      "</div>"
    );
  }

  function buttonLink(url, label) {
    return (
      '<a class="md-button md-button--primary" href="' +
      escapeAttr(url) +
      '">' +
      escapeHtml(label) +
      "</a>"
    );
  }

  function setRecommended(els, data) {
    if (!els.recWrap) {
      return;
    }
    var plat = detectPlatform();
    var osValue = plat.installerKey;
    var asset = data.installers && data.installers[osValue];

    if (asset) {
      els.recWrap.setAttribute("data-state", "ready");
      els.recTitle.textContent = "TASSEL 5 for " + plat.shortLabel;
      setRecMeta(els.recMeta, plat.note, asset.name, "memory");
      els.recBtn.setAttribute("href", asset.url);
      els.recBtn.textContent = "Download for " + plat.shortLabel;

      // Pre-select the detected OS in the picker.
      if (els.dist.value === "installer") {
        selectValue(els.os, osValue);
        els.os.dispatchEvent(new Event("change"));
      }
    } else {
      els.recWrap.setAttribute("data-state", "ready");
      els.recTitle.textContent = "TASSEL 5 Desktop";
      els.recMeta.textContent =
        "No native installer detected for your system \u2014 pick one below.";
      els.recBtn.setAttribute("href", "#choose");
      els.recBtn.textContent = "Choose a download";
    }
  }

  function setStandaloneRecommended(els, data) {
    if (!els.cliWrap) {
      return;
    }
    var latest = data.standaloneVersions[0];
    var entry = latest && data.standalone[latest];
    var primary = entry && (entry.targz || entry.zip);

    if (primary) {
      els.cliWrap.setAttribute("data-state", "ready");
      els.cliTitle.textContent = "TASSEL 5 Standalone v" + latest;
      setRecMeta(
        els.cliMeta,
        "Command line \u00b7 requires Java 21+",
        primary.name,
        "terminal"
      );
      els.cliBtn.setAttribute("href", primary.url);
      els.cliBtn.textContent = primary === entry.targz ? "Download .tar.gz" : "Download .zip";
    } else {
      els.cliWrap.setAttribute("data-state", "ready");
      els.cliTitle.textContent = "TASSEL 5 Standalone";
      els.cliMeta.textContent =
        "No standalone build available \u2014 pick a version below.";
      els.cliBtn.setAttribute("href", "#choose");
      els.cliBtn.textContent = "Choose a version";
    }
  }

  function detectPlatform() {
    var ua = navigator.userAgent || "";
    var platform = (navigator.platform || "").toLowerCase();

    if (/windows|win32|win64/i.test(ua) || platform.indexOf("win") === 0) {
      return {
        installerKey: "win-x64",
        label: "Windows (64-bit)",
        shortLabel: "Windows",
        note: ""
      };
    }

    if (/macintosh|mac os x/i.test(ua) || platform.indexOf("mac") === 0) {
      var arm = isAppleSilicon();
      return arm
        ? {
            installerKey: "mac-arm64",
            label: "macOS (Apple Silicon)",
            shortLabel: "macOS",
            note: "Apple Silicon"
          }
        : {
            installerKey: "mac-x64",
            label: "macOS (Intel)",
            shortLabel: "macOS",
            note: "Intel"
          };
    }

    return {
      installerKey: "linux-x64",
      label: "Linux / Unix",
      shortLabel: "Linux",
      note: ""
    };
  }

  // Best-effort Apple Silicon detection. Defaults to Apple Silicon on modern
  // macOS since it cannot be read directly from the browser.
  function isAppleSilicon() {
    try {
      var canvas = document.createElement("canvas");
      var gl =
        canvas.getContext("webgl") || canvas.getContext("experimental-webgl");
      if (gl) {
        var ext = gl.getExtension("WEBGL_debug_renderer_info");
        if (ext) {
          var renderer = gl.getParameter(ext.UNMASKED_RENDERER_WEBGL) || "";
          if (/intel/i.test(renderer)) {
            return false;
          }
          if (/apple/i.test(renderer)) {
            return true;
          }
        }
      }
    } catch (e) {
      /* ignore */
    }
    return true;
  }

  function setRecMeta(el, note, filename, noteIcon) {
    if (!el) {
      return;
    }
    if (filename) {
      var icon = noteIcon || "memory";
      var fileRow =
        '<span class="dl-hero__file-row">' +
        '<span class="material-symbols-outlined dl-hero__file-icon" aria-hidden="true">description</span>' +
        '<span class="dl-hero__file-viewport" title="' +
        escapeAttr(filename) +
        '">' +
        '<span class="dl-hero__file-track">' +
        '<span class="dl-hero__file-text">' +
        escapeHtml(filename) +
        "</span></span></span></span>";
      var noteRow = note
        ? '<span class="dl-hero__note-row">' +
          '<span class="material-symbols-outlined dl-hero__note-icon" aria-hidden="true">' +
          escapeHtml(icon) +
          "</span>" +
          '<span class="dl-hero__note">' +
          escapeHtml(note) +
          "</span></span>"
        : "";
      el.innerHTML = noteRow + fileRow;
      applyFileTicker(el);
      return;
    }
    el.textContent = note || "";
  }

  function applyFileTicker(metaEl) {
    var viewport = metaEl.querySelector(".dl-hero__file-viewport");
    if (!viewport) {
      return;
    }

    var track = viewport.querySelector(".dl-hero__file-track");
    var text = viewport.querySelector(".dl-hero__file-text");
    if (!track || !text) {
      return;
    }

    track.querySelectorAll('.dl-hero__file-text[aria-hidden="true"]').forEach(function (node) {
      node.remove();
    });
    viewport.classList.remove("dl-hero__file--ticker");
    track.style.removeProperty("--ticker-duration");

    if (text.scrollWidth <= viewport.clientWidth) {
      return;
    }

    var clone = text.cloneNode(true);
    clone.setAttribute("aria-hidden", "true");
    track.appendChild(clone);

    var duration = Math.max(8, Math.round(text.scrollWidth / 28));
    track.style.setProperty("--ticker-duration", duration + "s");
    viewport.classList.add("dl-hero__file--ticker");
  }

  function observeFileTicker(metaEl) {
    if (!metaEl || metaEl._tickerObserved || typeof ResizeObserver === "undefined") {
      return;
    }
    metaEl._tickerObserved = true;
    new ResizeObserver(function () {
      applyFileTicker(metaEl);
    }).observe(metaEl);
  }

  /* ---------- small DOM helpers ---------- */

  function fillSelect(select, options) {
    select.innerHTML = "";
    options.forEach(function (opt) {
      var o = document.createElement("option");
      o.value = opt.value;
      o.textContent = opt.label;
      select.appendChild(o);
    });
  }

  function selectValue(select, value) {
    for (var i = 0; i < select.options.length; i++) {
      if (select.options[i].value === value) {
        select.selectedIndex = i;
        return;
      }
    }
  }

  function labelFor(list, value) {
    for (var i = 0; i < list.length; i++) {
      if (list[i].value === value) {
        return list[i].label;
      }
    }
    return value;
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
