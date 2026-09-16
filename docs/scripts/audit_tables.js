// Paste into the browser console on any page of the served site (mkdocs serve).
// Loads every page of the nav in a hidden frame and lists table columns that
// hold prose but ended up narrower than MIN_PX. Run at the window width you
// care about; Material lays tables out by content, so a long unbreakable token
// in one column starves the others.
(async () => {
  const MIN_PX = 170;
  const links = [...document.querySelectorAll('.md-nav--primary a.md-nav__link[href]')]
    .map(a => a.href.split('#')[0]).filter((v, i, a) => a.indexOf(v) === i && v.startsWith(location.origin));
  const bad = []; let tables = 0;
  for (const url of links) {
    const f = document.createElement('iframe');
    f.style.cssText = `position:fixed;left:0;top:0;width:${innerWidth}px;height:900px;visibility:hidden`;
    document.body.appendChild(f);
    await new Promise(r => { f.onload = r; f.src = url; });
    await new Promise(r => setTimeout(r, 150));
    f.contentDocument.querySelectorAll('.md-typeset table:not([class])').forEach((t, ti) => {
      tables++;
      const row = t.tBodies[0] && t.tBodies[0].rows[0]; if (!row) return;
      const cells = [...row.cells], widths = cells.map(c => Math.round(c.getBoundingClientRect().width));
      cells.forEach((c, ci) => {
        const txt = c.textContent.trim();
        const code = [...c.querySelectorAll('code')].map(x => x.textContent).join('').trim();
        const prose = txt.length >= 40 && code.length < txt.length * 0.6;
        if (prose && widths[ci] < MIN_PX) bad.push(`${url.replace(location.origin, '')} table ${ti} column ${ci}: ${widths[ci]}px [${widths.join('/')}]`);
      });
    });
    f.remove();
  }
  console.log(`${links.length} pages, ${tables} tables, ${bad.length} narrow prose columns`);
  bad.forEach(b => console.log('  ' + b));
})();
