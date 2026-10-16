/* Command builder on "Choosing your inputs". The reader answers a few questions and gets a
   nextflow command. The rules follow the launch checks in main.nf, cited next to each.
   Page links come from hidden anchors in the page, so MkDocs resolves them for any URL style. */
(function () {
  "use strict";

  const DISCOVERY = [
    { key: "assemblies", flag: "--assemblies", sheet: "assemblies.csv", label: "Genome assemblies, one FASTA per haplotype" },
    { key: "pav", flag: "--pav", sheet: "pav.csv", label: "Phased assemblies to call with PAV" },
    { key: "longreads", flag: "--longreads", sheet: "longreads.csv", label: "Long reads, FASTQ or FASTA" },
    { key: "bams", flag: "--bams", sheet: "bams.csv", label: "Long reads aligned to the reference, BAM" },
    { key: "svs", flag: "--svs", sheet: "svs.csv", label: "Your own SV calls, one VCF per sample" },
  ];
  const VG = ["giraffe", "graphaligner"];

  const questions = (s) => [
    {
      key: "start", type: "radio", title: "What do you start from?",
      options: [
        { value: "discovery", label: "Assemblies, long reads, or SV calls per sample" },
        { value: "vcf", label: "One merged, sequence-resolved SV VCF" },
        { value: "rmdir", label: "The 2_Repeat_Filtering/ directory of an earlier run" },
        { value: "gvcf", label: "The pangenome.vcf of an earlier run" },
      ],
    },
    s.start === "discovery" && {
      key: "inputs", type: "check", title: "Which of these do you have? Tick all that apply.",
      options: DISCOVERY.map((d) => ({ value: d.key, label: d.label })),
    },
    {
      key: "human", type: "radio", title: "Do you want the human filter for young Alu, L1, SVA and HERV-K insertions?",
      options: [{ value: "no", label: "No" }, { value: "yes", label: "Yes, --human" }],
    },
    s.start !== "gvcf" && {
      key: "genotype", type: "radio", title: "Genotype read sets on the graph?",
      options: [{ value: "yes", label: "Yes" }, { value: "no", label: "No, stop after Stage 2" }],
    },
    genotyping(s) && {
      key: "method", type: "radio", title: "How should GraffiTE genotype them?",
      options: [
        { value: "pangenie", label: "PanGenie, the default: counts k-mers in short reads, no alignment" },
        { value: "giraffe", label: "vg giraffe, then vg call: short or long reads" },
        { value: "graphaligner", label: "GraphAligner, then vg call: long reads" },
      ],
    },
    genotyping(s) && VG.includes(s.method) && {
      key: "epi", type: "radio", title: "Map methylation from the genotyping BAMs (MM/ML tags)?",
      options: [{ value: "no", label: "No" }, { value: "yes", label: "Yes, --epigenomes" }],
    },
    {
      key: "profile", type: "radio", title: "Where will it run?",
      options: [
        { value: "standard", label: "This machine" },
        { value: "cluster", label: "A SLURM cluster" },
        { value: "cloud", label: "AWS Batch, with your own AWS settings" },
      ],
    },
  ].filter(Boolean);

  function genotyping(s) { return s.start === "gvcf" || s.genotype === "yes"; }

  // A note is a list of parts: plain strings, or [link text, page key, anchor].
  function build(s, rev) {
    const args = [], notes = [];
    const human = s.human === "yes";
    args.push(`-profile ${s.profile}`);
    if (s.profile === "cloud") {
      args.push("-c aws.config", "-w s3://your-bucket/work");
      notes.push(["aws.config holds your AWS Batch queue (process.queue) and region (aws.region). -profile cloud runs the image with Docker on AWS Batch; we have not tested it. ", ["Execution profiles", "install", "#execution-profiles"]]);
    }
    args.push("--reference reference.fa");
    // Stage 2 reads the TE library unless the run enters after RepeatMasker; hervk_annotate reads it under --human (main.nf:144,175).
    if (s.start !== "gvcf" && (s.start !== "rmdir" || human)) args.push("--TE_library TE_library.fa");
    if (s.start === "discovery") {
      const picked = DISCOVERY.filter((d) => s.inputs.includes(d.key));
      if (!picked.length) return { error: "Tick at least one input." };
      picked.forEach((d) => args.push(`${d.flag} ${d.sheet}`));
      const sheets = ["Samplesheet columns: "];
      picked.forEach((d, i) => sheets.push(i ? ", " : "", [d.flag, "sheets", "#" + d.flag.slice(1)]));
      notes.push(sheets.concat("."));
      if (s.inputs.includes("pav") && s.profile === "cloud") notes.push(["pav_asm cannot run under -profile cloud: its image is a Sylabs library image."]);
      if (s.inputs.includes("pav")) notes.push(["PAV 2 runs in its own container and is slow. We are working to add PAV 3, which is much faster. ", ["PAV", "discovery", "#from-assemblies-with-pav-pav"]]);
      if (s.inputs.includes("svs")) notes.push(["GraffiTE does not filter your VCFs: keep only insertions and deletions with sequence-resolved alleles. ", ["--svs", "skipping", "#-svs-your-own-per-sample-calls"]]);
    } else if (s.start === "vcf") {
      args.push("--vcf my_svs.vcf.gz");
      notes.push(["GraffiTE refuses a discovery flag beside --vcf. Split multi-allelic records first with bcftools norm -m-. ", ["--vcf", "skipping", "#-vcf-one-merged-vcf"]]);
    } else if (s.start === "rmdir") {
      args.push("--RM_dir earlier/out/2_Repeat_Filtering");
      notes.push(["The run restarts at the TSD search, so new --tsd_win, trusted, --human or HERV-K settings take effect. ", ["--RM_dir", "skipping", "#-rm_dir-repeatmasker-output-of-an-earlier-run"]]);
    } else {
      args.push("--graffite_vcf earlier/out/3_TSD_search/pangenome.vcf");
      notes.push(["You can reuse the earlier graph, alignments or vg call VCFs too. ", ["Resuming and skipping work", "skipping", "#-graffite_vcf-a-finished-pangenomevcf"]]);
    }
    if (human) args.push("--human");
    if (!genotyping(s)) {
      args.push("--genotype false");
      return { cmd: command(rev, args), notes };
    }
    args.push("--genotype_with reads.csv", `--graph_method ${s.method}`);
    notes.push(["reads.csv lists sample, path and type, one read set per row: FASTQ or BAM, short or long reads. ", ["--genotype_with", "sheets", "#-genotype_with"]]);
    // Two launch checks stop a --human run here: main.nf:58-63 and main.nf:65-71.
    if (human && s.start === "gvcf") {
      args.push("--hervk_reconcile false");
      notes.push(["With --graffite_vcf, GraffiTE skips the step that writes the HERV-K tables and will not start unless --hervk_reconcile false is set, so HERV-K loci are not consolidated. Start from --RM_dir to keep that step."]);
    } else if (human && s.method === "pangenie") {
      args.push("--hervk_reconcile false");
      notes.push(["GraffiTE cannot consolidate PanGenie genotypes per HERV-K locus and will not start --human with PanGenie unless --hervk_reconcile false is set. Choose vg giraffe to keep that step. ", ["Human MEIs", "human"]]);
    }
    if (VG.includes(s.method) && s.epi === "yes") {
      args.push("--epigenomes");
      notes.push(["Every row of reads.csv must be a BAM with MM/ML tags: GraffiTE genotypes a FASTQ row, then leaves it out of the merged VCF without a warning. ", ["Epigenetic mapping", "methylation"]]);
    }
    return { cmd: command(rev, args), notes };
  }

  function command(rev, args) {
    return `nextflow run cgroza/GraffiTE -r ${rev} -latest \\\n  ` + args.join(" \\\n  ");
  }

  // Flags go in <code>: in the body font, "--" reads as a single dash.
  function rich(text) {
    const f = document.createDocumentFragment();
    text.split(/(--[A-Za-z_]+)/).forEach((part, i) => part && f.append(i % 2 ? el("code", { text: part }) : part));
    return f;
  }

  function el(tag, attrs, ...kids) {
    const e = document.createElement(tag);
    Object.entries(attrs || {}).forEach(([k, v]) => (k === "text" ? (e.textContent = v) : e.setAttribute(k, v)));
    kids.forEach((k) => k && e.append(k));
    return e;
  }

  function init(root) {
    if (root.dataset.ready) return;
    root.dataset.ready = "1";
    const rev = root.dataset.revision;
    const links = {};
    root.querySelectorAll("a[data-key]").forEach((a) => (links[a.dataset.key] = a.getAttribute("href")));
    const s = { start: "discovery", inputs: ["assemblies"], human: "no", genotype: "yes", method: "pangenie", epi: "no", profile: "standard" };
    const form = el("form", { class: "gt-picker-form" });
    const out = el("div", { class: "gt-picker-out", "aria-live": "polite" });
    form.addEventListener("submit", (ev) => ev.preventDefault());
    root.replaceChildren(form, out);

    function render() {
      form.replaceChildren();
      questions(s).forEach((q, qi) => {
        const fs = el("fieldset", {}, el("legend", {}, rich(`${qi + 1}. ${q.title}`)));
        q.options.forEach((o) => {
          const id = `gt-${q.key}-${o.value}`;
          const input = el("input", { type: q.type === "check" ? "checkbox" : "radio", name: q.key, id, value: o.value });
          input.checked = q.type === "check" ? s[q.key].includes(o.value) : s[q.key] === o.value;
          input.addEventListener("change", () => {
            s[q.key] = q.type === "check"
              ? [...form.querySelectorAll(`input[name=${q.key}]:checked`)].map((i) => i.value)
              : o.value;
            render();
            document.getElementById(id)?.focus();
          });
          fs.append(el("label", { for: id }, input, el("span", {}, rich(o.label))));
        });
        form.append(fs);
      });
      const r = build(s, rev);
      out.replaceChildren();
      if (r.error) { out.append(el("p", { class: "gt-picker-error", text: r.error })); return; }
      const code = el("code", { text: r.cmd });
      const copy = el("button", { type: "button", class: "md-button gt-picker-copy", text: "Copy" });
      copy.addEventListener("click", () => {
        const select = () => {
          const range = document.createRange(); range.selectNodeContents(code);
          const sel = getSelection(); sel.removeAllRanges(); sel.addRange(range);
        };
        const done = () => { copy.textContent = "Copied"; setTimeout(() => (copy.textContent = "Copy"), 1500); };
        if (navigator.clipboard) navigator.clipboard.writeText(r.cmd).then(done, select); else select();
      });
      out.append(el("pre", {}, code), copy);
      if (r.notes.length) {
        const ul = el("ul", {});
        r.notes.forEach((parts) => {
          const li = el("li", {});
          parts.forEach((p) => {
            if (typeof p === "string") li.append(rich(p));
            else if (links[p[1]]) li.append(el("a", { href: links[p[1]] + (p[2] || "") }, rich(p[0])));
            else li.append(rich(p[0]));
          });
          ul.append(li);
        });
        out.append(ul);
      }
    }
    render();
  }

  const run = () => document.querySelectorAll(".gt-picker").forEach(init);
  // Material's instant navigation swaps pages without a reload and emits each one on document$.
  if (typeof document$ !== "undefined") document$.subscribe(run);
  else document.addEventListener("DOMContentLoaded", run);
})();
