/*
 * Shared GitHub Releases feed for the TASSEL docs site.
 *
 * The download page and the nightly builds page both need the full release
 * list. Unauthenticated GitHub API calls are limited to 60 per hour per IP, so
 * the response is fetched once, cached in sessionStorage, and shared through
 * window.tasselReleases.
 */
(function () {
  "use strict";

  var REPO = "maize-genetics/tassel";
  var API = "https://api.github.com/repos/" + REPO + "/releases?per_page=100";
  var CACHE_KEY = "tassel-releases-v1";
  var CACHE_TTL_MS = 30 * 60 * 1000; // 30 minutes

  // In-flight promise, so two pages initialising in the same tick issue one
  // request rather than two.
  var pending = null;

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

    if (pending) {
      return pending;
    }

    pending = fetch(API, { headers: { Accept: "application/vnd.github+json" } })
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
        pending = null;
        return releases;
      })
      .catch(function (err) {
        pending = null;
        throw err;
      });

    return pending;
  }

  window.tasselReleases = {
    repo: REPO,
    releasesUrl: "https://github.com/" + REPO + "/releases",
    fetch: fetchReleases
  };
})();
