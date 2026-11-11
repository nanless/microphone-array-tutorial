(() => {
  const scene = document.querySelector(".course-demo");
  const switches = scene.querySelector(".course-switch");
  const arrival = scene.querySelector(".course-arrival");
  const sources = [...scene.querySelectorAll("[data-course-source]")];
  const reduced = window.matchMedia("(prefers-reduced-motion: reduce)");
  const replay = scene.querySelector("[data-course-replay]");
  switches.hidden = false;

  function animate() {
    scene.classList.remove("course-playing");
    if (reduced.matches) return;
    // Resolve the removed state before restarting the finite animation.
    void scene.offsetWidth;
    scene.classList.add("course-playing");
  }

  function motionPreferenceChanged() {
    replay.hidden = reduced.matches;
    animate();
  }

  for (const button of sources) {
    button.addEventListener("click", () => {
      const side = button.dataset.courseSource;
      scene.dataset.source = side;
      for (const option of sources) option.setAttribute("aria-pressed", String(button === option));
      arrival.textContent = side === "left" ? "左边 A 发声：1 → 2 → 3 → 4" : "右边 B 发声：4 → 3 → 2 → 1";
      animate();
    });
  }
  replay.addEventListener("click", animate);
  reduced.addEventListener("change", motionPreferenceChanged);
  motionPreferenceChanged();

  const clips = [...document.querySelectorAll(".course-audio-grid audio")];
  for (const clip of clips) {
    clip.addEventListener("play", () => {
      for (const other of clips) if (other !== clip) other.pause();
    });
  }
})();
