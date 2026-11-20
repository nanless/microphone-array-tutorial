(() => {
  for (const peer of document.querySelectorAll("a[data-lang-peer]")) {
    if (location.protocol === "file:" || location.hostname === "localhost" || location.hostname === "127.0.0.1") {
      peer.href = peer.dataset.langPeer;
    }
  }
  const scene = document.querySelector(".course-demo");
  const switches = scene.querySelector(".course-switch");
  const arrival = scene.querySelector(".course-arrival");
  const sources = [...scene.querySelectorAll("[data-course-source]")];
  const reduced = window.matchMedia("(prefers-reduced-motion: reduce)");
  const replay = scene.querySelector("[data-course-replay]");
  const modeNotes = {
    mvdr: "Active mode: MVDR locks onto −35° target A (0 dB distortionless) while carving a deep spatial null toward +30° interferer B.",
    left: "Active mode: Delay-and-Sum (DSB) compensates A's arrival delays (M1→M4) with equal weights (WNG = 6.02 dB).",
    right: "Active mode: Re-steered toward +30° source B (arrival order flips M4→M3→M2→M1); B sums coherently while A cancels.",
  };
  switches.hidden = false;

  function animate() {
    scene.classList.remove("course-playing");
    if (reduced.matches) return;
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
      arrival.textContent = modeNotes[side] || modeNotes.mvdr;
      animate();
    });
  }
  replay.addEventListener("click", animate);
  reduced.addEventListener("change", motionPreferenceChanged);
  motionPreferenceChanged();

  const localization = document.querySelector(".course-localization");
  const scanSwitch = localization.querySelector(".course-scan-switch");
  const scanResult = localization.querySelector(".course-scan-result");
  const scanOptions = [...scanSwitch.querySelectorAll("[data-course-scan]")];
  const scanLabels = {
    a: "A · Est. 35°, Optical Truth 34.4°",
    b: "B · Est. 16°, Optical Truth 18.5°",
  };
  scanSwitch.hidden = false;
  for (const button of scanOptions) {
    button.addEventListener("click", () => {
      const source = button.dataset.courseScan;
      localization.dataset.scan = source;
      for (const option of scanOptions) option.setAttribute("aria-pressed", String(option === button));
      scanResult.textContent = scanLabels[source];
    });
  }

  const clips = [...document.querySelectorAll(".course-audio-grid audio")];
  for (const clip of clips) {
    clip.addEventListener("play", () => {
      for (const other of clips) if (other !== clip) other.pause();
    });
  }

  const jumps = document.querySelector(".course-audio-jumps");
  const seekStatus = jumps.querySelector(".course-seek-status");
  const pending = new Set();
  let seekSeconds;

  function seekWhenReady(clip) {
    if (!pending.has(clip)) return;
    const target = Number(seekSeconds);
    const ranges = clip.seekable;
    for (let i = 0; i < ranges.length; i += 1) {
      if (ranges.start(i) <= target && target < ranges.end(i)) {
        clip.currentTime = target;
        pending.delete(clip);
        if (pending.size === 0) seekStatus.textContent = `Synchronized all players to ${seekSeconds} s — press Play on any track to compare.`;
        return;
      }
    }
  }

  for (const clip of clips) {
    for (const event of ["loadedmetadata", "progress", "canplay", "canplaythrough"]) {
      clip.addEventListener(event, () => seekWhenReady(clip));
    }
  }

  jumps.hidden = false;
  for (const button of jumps.querySelectorAll("[data-course-seek]")) {
    button.addEventListener("click", () => {
      seekSeconds = button.dataset.courseSeek;
      seekStatus.textContent = `Seeking to ${seekSeconds} s…`;
      for (const clip of clips) pending.add(clip);
      for (const clip of clips) {
        clip.pause();
        if (clip.readyState === 0 && clip.networkState !== HTMLMediaElement.NETWORK_LOADING) {
          clip.preload = "auto";
          clip.load();
        }
        seekWhenReady(clip);
      }
    });
  }
})();
