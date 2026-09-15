#!/usr/bin/env bash
# Regression test for truvari_merge with more than one input VCF.
#
# The sharded collapse (e477ec1) wrote truvari_merged.vcf.gz, but the next
# line still read truvari_merged.vcf, so a run that reached truvari_merge
# with two or more VCFs failed there with "Failed to open file". The test
# takes the script straight from module/main.nf, runs it on two svim-asm
# style VCFs and checks SVs.vcf.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
here=$PWD
module="$here/../../module/main.nf"
export PATH="$(cd ../../bin && pwd):$PATH"

for t in bcftools bgzip tabix truvari; do
  command -v $t >/dev/null || { echo "  [skip] $t not on PATH"; exit 0; }
done
python3 -c 'import vcfpy' 2>/dev/null || { echo "  [skip] vcfpy not importable"; exit 0; }

fail=0
chk(){ if [[ "$2" == "$3" ]]; then echo "  [ ok ] $1"; else
       echo "  [FAIL] $1: got '$2', want '$3'"; fail=1; fi; }

tmp=$(mktemp -d); trap 'rm -rf "$tmp"' EXIT
cd "$tmp"

# shorten_ids.py names /usr/bin/python3, which may not be the python3 that
# has vcfpy. The script calls it by name, so a shim first on PATH redirects it.
mkdir shim
printf '#!/usr/bin/env bash\nexec python3 %s "$@"\n' "$(command -v shorten_ids.py)" > shim/shorten_ids.py
chmod +x shim/shorten_ids.py
export PATH="$tmp/shim:$PATH"

python3 - "$module" <<'EOF'
import random, re, sys
random.seed(3)
seq = ''.join(random.choice('ACGT') for _ in range(20000))
open('ref.fa', 'w').write('>chr1\n' + '\n'.join(seq[i:i+60] for i in range(0, len(seq), 60)) + '\n')
def s(n): return ''.join(random.choice('ACGT') for _ in range(n))
shared, a_only, b_only = s(300), s(250), s(400)
def hdr(name):
    return ['##fileformat=VCFv4.2', '##contig=<ID=chr1,length=20000>',
            '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="type">',
            '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="length">',
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
            '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t' + name]
def ins(p, i, x): return f"chr1\t{p}\t{i}\t{seq[p-1]}\t{seq[p-1]+x}\t.\tPASS\tSVTYPE=INS;SVLEN={len(x)}\tGT\t1/1"
def dl(p, i, n): return f"chr1\t{p}\t{i}\t{seq[p-1:p+n]}\t{seq[p-1]}\t.\tPASS\tSVTYPE=DEL;SVLEN=-{n}\tGT\t1/1"
# The two insertions near 1000 differ by 2 bp of position and 3 bp of length,
# so truvari collapses them; the deletions are identical.
open('A.vcf', 'w').write('\n'.join(hdr('A') + [ins(1000, 'A.INS.1', shared), ins(5000, 'A.INS.2', a_only), dl(9000, 'A.DEL.1', 500)]) + '\n')
open('B.vcf', 'w').write('\n'.join(hdr('B') + [ins(1002, 'B.INS.1', shared[:-3]), ins(15000, 'B.INS.2', b_only), dl(9000, 'B.DEL.1', 500)]) + '\n')

# Nextflow's interpolation, done by hand for this process.
body = re.search(r'process truvari_merge \{.*?script:\n  """\n(.*?)\n  """', open(sys.argv[1]).read(), re.S).group(1)
body = (body.replace('\\$', '$').replace('${vcfs}', 'A.vcf.gz B.vcf.gz')
            .replace('${task.cpus}', '2').replace('${ref}', 'ref.fa').replace('${from_vcf}', 'false'))
open('truvari_merge.sh', 'w').write(body)
EOF
for s in A B; do bcftools sort -Oz -o $s.vcf.gz $s.vcf 2>/dev/null; done

# BSD xargs gives -I replacements 255 bytes by default, and this command is
# longer. Where xargs is not GNU (macOS), the shards run one after the other.
if ! xargs --version >/dev/null 2>&1; then
  python3 - <<'EOF'
b = open('truvari_merge.sh').read()
b = b.replace('xargs -P "2" -I{} bash -c \'', 'while read -r s; do bash -c \'').replace('shard="{}"', 'shard="$1"')
b = b.replace('rm "collapsed/${base}.vcf"\n    \'\n', 'rm "collapsed/${base}.vcf"\n    \' _ "$s"; done\n')
open('truvari_merge.sh', 'w').write(b)
EOF
fi

rc=0; bash -ue truvari_merge.sh > run.log 2>&1 || rc=$?
chk "script exits 0" "$rc" "0"
[[ -f SVs.vcf ]] || { echo "  [FAIL] no SVs.vcf"; tail -5 run.log; exit 1; }
chk "records after collapse"   "$(bcftools view -H SVs.vcf | wc -l | tr -d ' ')" "4"
chk "collapsed insertion, both genotypes" "$(bcftools view -H SVs.vcf | awk '$2==1000' | cut -f10,11)" "1/1	1/1"
chk "missing GT set to 0/0"    "$(bcftools view -H SVs.vcf | awk '$2==5000' | cut -f10,11)" "1/1	0/0"
chk "SVLEN filled"             "$(bcftools view -H SVs.vcf | awk '$2==15000' | grep -o 'SVLEN=[-0-9]*')" "SVLEN=400"
chk "IDs shortened and unique" "$(bcftools view -H SVs.vcf | cut -f3 | sort -u | wc -l | tr -d ' ')" "4"

exit $fail
