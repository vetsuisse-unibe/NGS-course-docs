// NGS course exercises, HTML version 2: copy buttons on command blocks and
// highlighting of the current section in the table of contents.
document.addEventListener("DOMContentLoaded", function () {
  // Copy buttons on every command block
  document.querySelectorAll("pre.shell, pre.sourceCode.shell").forEach(function (pre) {
    var wrap = document.createElement("div");
    wrap.className = "terminal";
    pre.parentNode.insertBefore(wrap, pre);
    wrap.appendChild(pre);

    var btn = document.createElement("button");
    btn.className = "copy-btn";
    btn.type = "button";
    btn.textContent = "Copy";
    btn.addEventListener("click", function () {
      var text = pre.innerText.replace(/\n$/, "");
      var done = function () {
        btn.textContent = "Copied";
        btn.classList.add("copied");
        setTimeout(function () { btn.textContent = "Copy"; btn.classList.remove("copied"); }, 1500);
      };
      if (navigator.clipboard && window.isSecureContext) {
        navigator.clipboard.writeText(text).then(done, function () {});
      } else {
        var ta = document.createElement("textarea");
        ta.value = text;
        document.body.appendChild(ta);
        ta.select();
        try { document.execCommand("copy"); done(); } catch (e) {}
        document.body.removeChild(ta);
      }
    });
    wrap.appendChild(btn);
  });

  // Highlight the section being read in the table of contents
  var links = {};
  document.querySelectorAll("nav#TOC a").forEach(function (a) {
    links[decodeURIComponent(a.getAttribute("href").slice(1))] = a;
  });
  var sections = document.querySelectorAll("h2[id]");
  if (!("IntersectionObserver" in window) || sections.length === 0) return;
  var current = null;
  var observer = new IntersectionObserver(function (entries) {
    entries.forEach(function (entry) {
      if (entry.isIntersecting && links[entry.target.id]) {
        if (current) current.classList.remove("active");
        current = links[entry.target.id];
        current.classList.add("active");
      }
    });
  }, { rootMargin: "0px 0px -70% 0px" });
  sections.forEach(function (s) { observer.observe(s); });
});
