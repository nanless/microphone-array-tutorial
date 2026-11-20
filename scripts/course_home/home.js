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
    mvdr: "当前模式：MVDR 锁住 −35° 目标 A（0 dB 无失真），同时在 +30° 干扰 B 方向刻出空间零陷。",
    left: "当前模式：延迟求和（DSB）补齐 A 的到达时差（M1→M4），四路同相平均，白噪声增益 6.02 dB。",
    right: "当前模式：改向对准 +30° 声源 B（到达次序反转为 M4→M3→M2→M1），B 同相增强而 A 错位削弱。",
  };
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
  const scanLabels = {a: "A · 估计 35°，实测 34.4°", b: "B · 估计 16°，实测 18.5°"};
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
        if (pending.size === 0) seekStatus.textContent = `已同步定位到 ${seekSeconds} 秒，可点击播放任一路对比。`;
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
      seekStatus.textContent = `准备跳到 ${seekSeconds} 秒…`;
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
