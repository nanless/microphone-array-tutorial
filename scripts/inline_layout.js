/* Shared publication layout. TeX, code, links and attributes remain intact. */
(function () {
  "use strict";
  const mathClass = "tutorial-math-tail";
  const idClass = "tutorial-exercise-id";
  const ignored = "pre,code,a,script,style,textarea,mjx-container";
  const suffix = /^(\s*(?:(?:m\/s[²³]?|rad\/s|rad\/m|deg\/s|μs|µs|ns|ms|s|dBFS|dB|kHz|MHz|Hz|ppm|rad|cm|mm|m|V|deg|°|%|％)(?:[²³])?(?=$|[\s，。；、：！？）,.!?;:)\]\u4e00-\u9fff])\s*[，。；、：！？）,.;:!?\]]*|[，。；、：！？）,.;:!?\]]+))/;
  let ready = false;
  let scheduled = false;
  let printing = false;

  function roots() {
    const candidates = Array.from(document.querySelectorAll(".main,.chap"));
    return candidates.filter(root => !candidates.some(other => other !== root && other.contains(root)));
  }

  function unwrapMath(root) {
    for (const wrapper of root.querySelectorAll("." + mathClass)) {
      const children = Array.from(wrapper.childNodes);
      wrapper.replaceWith(...children);
      // Merge only our adjacent ordinary text, never normalize a MathJax subtree.
      for (const node of children) {
        if (node.nodeType !== Node.TEXT_NODE) continue;
        while (node.nextSibling && node.nextSibling.nodeType === Node.TEXT_NODE) {
          node.data += node.nextSibling.data;
          node.nextSibling.remove();
        }
      }
    }
  }

  function protectExerciseIds(root) {
    let count = 0;
    {
      const walker = document.createTreeWalker(root, NodeFilter.SHOW_TEXT);
      const nodes = [];
      let node;
      while ((node = walker.nextNode())) {
        if (!node.parentElement.closest(ignored + ",." + idClass)) nodes.push(node);
      }
      for (const text of nodes) {
        // A filename/path or a longer identifier is not a stable exercise ID.
        const pattern = /(^|[^\p{L}\p{N}\p{M}_/\\.-])(E\d{2}-\d{2})(?![\p{L}\p{N}\p{M}_/\\-]|\.[\p{L}\p{N}\p{M}])/gu;
        const matches = Array.from(text.data.matchAll(pattern));
        for (const match of matches.reverse()) {
          const start = match.index + match[1].length;
          text.splitText(start + match[2].length);
          const middle = text.splitText(start);
          const span = document.createElement("span");
          span.className = idClass;
          middle.before(span);
          span.append(middle);
          count++;
        }
      }
    }
    return count;
  }

  function markShortLeads(root) {
    let count = 0;
    for (const paragraph of root.querySelectorAll("p")) {
      if (paragraph.closest("pre,code,a,script,style,textarea")) continue;
      const significant = Array.from(paragraph.childNodes).filter(node =>
        node.nodeType !== Node.TEXT_NODE || node.data.trim());
      const next = paragraph.nextElementSibling;
      const isEmptyAnchor = node =>
        node.nodeType === Node.ELEMENT_NODE &&
        node.tagName === "A" && node.id &&
        !node.hasAttribute("href") &&
        Array.from(node.childNodes).every(child =>
          child.nodeType === Node.TEXT_NODE && !child.data.trim());
      const isAnchorAlias = node =>
        node.nodeType === Node.ELEMENT_NODE &&
        node.tagName === "SPAN" && node.id &&
        node.classList.contains("anchor-alias") &&
        Array.from(node.childNodes).every(child =>
          child.nodeType === Node.TEXT_NODE && !child.data.trim());
      const anchorOnly = significant.length >= 1 &&
        significant.some(isEmptyAnchor) &&
        significant.every(node => isEmptyAnchor(node) || isAnchorAlias(node));
      paragraph.classList.toggle("tutorial-anchor-only", !!anchorOnly);
      const labelOnly = significant.length === 1 || (
        significant.length === 2 &&
        significant[1].nodeType === Node.TEXT_NODE &&
        /^[：:.。]$/.test(significant[1].data.trim()) &&
        significant[0].textContent.trim().length > 0);
      const short = labelOnly &&
        significant[0].nodeType === Node.ELEMENT_NODE &&
        significant[0].tagName === "STRONG" &&
        paragraph.textContent.trim().length <= 60 &&
        next && /^(P|UL|OL|BLOCKQUOTE)$/.test(next.tagName) &&
        next.textContent.trim().length > 0;
      paragraph.classList.toggle("tutorial-short-lead", !!short);
      if (short) count++;

      const firstSig = significant[0];
      const firstStrong = firstSig && firstSig.nodeType === Node.ELEMENT_NODE &&
        firstSig.tagName === "STRONG" ? firstSig : null;
      const strongText = firstStrong ? firstStrong.textContent.trim() : "";
      const isExLead = !!(firstStrong && firstStrong.querySelector("." + idClass) &&
        /^E\d{2}-\d{2}\b/.test(strongText));
      const isSolutionLead = !isExLead &&
        /^(参考答案|解答|答案|解(?=[：:。.\s]|$))/.test(strongText);
      const isProblemLead = !isExLead && !isSolutionLead &&
        /^(题目|问题|思考题)/.test(strongText);
      paragraph.classList.toggle("tutorial-exercise-lead", isExLead);
      paragraph.classList.toggle("tutorial-solution-lead", isSolutionLead);
      paragraph.classList.toggle("tutorial-problem-lead", isProblemLead);
    }
    for (const heading of root.querySelectorAll("h3,h4")) {
      const hasExId = !!heading.querySelector("." + idClass);
      const text = heading.textContent.trim();
      const isWorked = !hasExId && heading.tagName === "H3" &&
        (/^(算例\s*\d|例\s*\d)/.test(text) || /完整手算/.test(text));
      heading.classList.toggle("tutorial-exercise-heading", hasExId);
      heading.classList.toggle("tutorial-worked-heading", isWorked);
    }
    for (const item of root.querySelectorAll("ol > li")) {
      const existingCard = item.querySelector(":scope > .tutorial-exercise-card");
      const scopeRoot = existingCard || item;
      const hasSolution = !!scopeRoot.querySelector(":scope > p.tutorial-solution-lead");
      const startsQuestion = /^思考题[：:]/.test(scopeRoot.textContent.trim());
      if ((hasSolution || startsQuestion) && !existingCard) {
        const card = document.createElement("article");
        card.className = "tutorial-exercise-card";
        card.append(...Array.from(item.childNodes));
        item.append(card);
      }
      item.classList.toggle("tutorial-exercise-item", hasSolution || startsQuestion);
    }
    return count;
  }

  function apply() {
    const content = roots();
    for (const root of content) unwrapMath(root);
    const plans = [];
    const hosts = new Map();
    let skippedWide = 0;
    let idCount = 0;
    let leadCount = 0;
    for (const root of content) {
      idCount += protectExerciseIds(root);
      leadCount += markShortLeads(root);
      // Read every width before any math wrapper is inserted, to avoid repeated
      // synchronous layout across thousands of formulas in the combined book.
      for (const math of root.querySelectorAll('mjx-container:not([display="true"])')) {
        if (math.closest("pre,code,a,script,style,textarea")) continue;
        const text = math.nextSibling;
        if (!text || text.nodeType !== Node.TEXT_NODE) continue;
        const match = text.data.match(suffix);
        if (!match || match[0].length > 16) continue;
        const host = math.closest("p,li,td,th,h1,h2,h3,h4");
        if (!host) continue;
        if (!hosts.has(host)) {
          const style = getComputedStyle(host);
          hosts.set(host, host.clientWidth - parseFloat(style.paddingLeft || 0) -
            parseFloat(style.paddingRight || 0));
        }
        const range = document.createRange();
        range.setStart(text, 0);
        range.setEnd(text, match[0].length);
        const style = getComputedStyle(math);
        const width = math.getBoundingClientRect().width + range.getBoundingClientRect().width +
          parseFloat(style.marginLeft || 0) + parseFloat(style.marginRight || 0);
        if (width > hosts.get(host) - 8) {
          skippedWide++;
          continue;
        }
        plans.push({math, text, length: match[0].length});
      }
    }
    for (const {math, text, length} of plans) {
      const wrapper = document.createElement("span");
      wrapper.className = mathClass;
      math.before(wrapper);
      wrapper.append(math);
      text.splitText(length);
      wrapper.append(text);
    }
    ready = true;
    const result = {grouped: plans.length, skippedWide, newExerciseIds: idCount, shortLeads: leadCount};
    window.ArrayTutorialLayout.lastRun = result;
    return result;
  }

  function schedule() {
    if (!ready || scheduled || printing) return;
    scheduled = true;
    requestAnimationFrame(() => {
      scheduled = false;
      if (printing) return;
      apply();
    });
  }

  window.ArrayTutorialLayout = {apply, lastRun: null};
  window.addEventListener("resize", schedule);
  // Chrome fires beforeprint while screen styles still apply. Remove screen
  // groups before they can enlarge the print layout. The combined book is
  // measured at its fixed A4 geometry during pageReady instead. Do not mutate
  // thousands of nodes after the browser has begun its print snapshot: this
  // can truncate the final pages even when a small fixture prints correctly.
  // A synthetic beforeprint event is not a print test.
  window.addEventListener("beforeprint", () => {
    printing = true;
    if (ready && document.documentElement.dataset.tutorialPrintLayout !== "a4") {
      for (const root of roots()) unwrapMath(root);
    }
  });
  window.addEventListener("afterprint", () => {
    printing = false;
    if (document.documentElement.dataset.tutorialPrintLayout !== "a4") schedule();
  });
  if (document.fonts) document.fonts.addEventListener("loadingdone", schedule);
})();
