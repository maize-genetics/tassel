(function () {
  "use strict";

  function initGallery() {
    var gallery = document.getElementById("gallery");
    var overlay = document.getElementById("gallery-lightbox");
    if (!gallery || !overlay) {
      return;
    }

    var triggers = Array.prototype.slice.call(
      gallery.querySelectorAll(".gallery-trigger")
    );
    if (triggers.length === 0) {
      return;
    }

    var items = triggers.map(function (trigger) {
      var img = trigger.querySelector("img");
      return {
        src: img ? img.getAttribute("src") : "",
        alt: img ? img.getAttribute("alt") || "" : "",
        caption: trigger.getAttribute("data-caption") || ""
      };
    });

    var overlayImg = overlay.querySelector(".lightbox-img");
    var overlayCaption = overlay.querySelector(".lightbox-caption");
    var overlayCounter = overlay.querySelector(".lightbox-counter");
    var closeBtn = overlay.querySelector(".lightbox-close");
    var prevBtn = overlay.querySelector(".lightbox-prev");
    var nextBtn = overlay.querySelector(".lightbox-next");

    var currentIndex = 0;

    function render() {
      var item = items[currentIndex];
      overlayImg.setAttribute("src", item.src);
      overlayImg.setAttribute("alt", item.alt);
      overlayCaption.textContent = item.caption;
      overlayCounter.textContent = currentIndex + 1 + " / " + items.length;
    }

    function open(index) {
      currentIndex = index;
      render();
      overlay.classList.add("open");
      overlay.setAttribute("aria-hidden", "false");
      document.body.classList.add("lightbox-open");
      closeBtn.focus();
    }

    function close() {
      overlay.classList.remove("open");
      overlay.setAttribute("aria-hidden", "true");
      document.body.classList.remove("lightbox-open");
    }

    function show(index) {
      currentIndex = (index + items.length) % items.length;
      render();
    }

    triggers.forEach(function (trigger, index) {
      trigger.addEventListener("click", function () {
        open(index);
      });
    });

    prevBtn.addEventListener("click", function () {
      show(currentIndex - 1);
    });

    nextBtn.addEventListener("click", function () {
      show(currentIndex + 1);
    });

    closeBtn.addEventListener("click", close);

    overlay.addEventListener("click", function (event) {
      if (event.target === overlay) {
        close();
      }
    });

    document.addEventListener("keydown", function (event) {
      if (!overlay.classList.contains("open")) {
        return;
      }
      if (event.key === "Escape") {
        close();
      } else if (event.key === "ArrowLeft") {
        show(currentIndex - 1);
      } else if (event.key === "ArrowRight") {
        show(currentIndex + 1);
      }
    });
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", initGallery);
  } else {
    initGallery();
  }
})();
