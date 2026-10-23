"""Build the flat SNAPP upload for The Mathematical Intelligencer from main.tex:
no subfolders, figures renamed Fig<n><part>.pdf in order of appearance, source + compiled PDF + zip.
usage: python3.12 make_submission.py   -> submission/ and submission.zip"""
import os, re, shutil, subprocess, zipfile
root = os.path.dirname(os.path.abspath(__file__)); out = os.path.join(root, "submission")
FIGS = {  # source figure -> submission name, in the order the figures appear
    "parity_k150.pdf": "Fig1a.pdf", "parity_k1000.pdf": "Fig1b.pdf",
    "parity_res_N256_t003.pdf": "Fig2.pdf", "parity_late_N256.pdf": "Fig3b.pdf",
    **{f"glyco_{k}_t{t:04d}.pdf": f"Fig4{p}_{k}.pdf" for p, t in zip("abcd", (0, 200, 350, 1000)) for k in ("field", "parity")},
}
SHARED = {"parity_k2000.pdf": ("Fig1c.pdf", "Fig3a.pdf")}  # one source file used in two figures

tex = open(os.path.join(root, "main.tex")).read()
tex = tex.replace("\\graphicspath{{figures/}}\n", "")
for src, dst in FIGS.items():
    assert tex.count("{" + src + "}") == 1, src
    tex = tex.replace("{" + src + "}", "{" + dst + "}")
src = "{parity_k2000.pdf}"; a, b = SHARED["parity_k2000.pdf"]
assert tex.count(src) == 2
tex = tex.replace(src, "{" + a + "}", 1).replace(src, "{" + b + "}", 1)
left = re.findall(r"\\includegraphics(?:\[[^\]]*\])?\{([^}]+)\}", tex)
assert all(f.startswith("Fig") for f in left), left

shutil.rmtree(out, ignore_errors=True); os.makedirs(out)
open(os.path.join(out, "main.tex"), "w").write(tex)
for s, d in FIGS.items(): shutil.copy(os.path.join(root, "figures", s), os.path.join(out, d))
for d in SHARED["parity_k2000.pdf"]: shutil.copy(os.path.join(root, "figures", "parity_k2000.pdf"), os.path.join(out, d))
for _ in range(2):
    r = subprocess.run(["pdflatex", "-interaction=nonstopmode", "main.tex"], cwd=out, capture_output=True, text=True)
assert r.returncode == 0 and "undefined" not in r.stdout, r.stdout[-2000:]
for ext in ("aux", "log", "out"): os.remove(os.path.join(out, "main." + ext))
with zipfile.ZipFile(os.path.join(root, "submission.zip"), "w", zipfile.ZIP_DEFLATED) as z:
    for f in sorted(os.listdir(out)): z.write(os.path.join(out, f), f)
print("\n".join(sorted(os.listdir(out))))
