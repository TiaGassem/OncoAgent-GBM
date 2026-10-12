"""Labelled interactive 3D docking viewer (3Dmol.js, pinned version) for Streamlit.

Shows: protein cartoon + name, ligand pose(s) with Vina energy, pose selector, pocket residues
(name+number), dashed distance lines with Å labels for N/O pairs <= 3.5 Å (candidates only),
grid box, optional crystal reference ligand (for redocking overlay). Spin is OFF by default.
All geometry is computed in Python from the actual files; nothing is hard-coded.
"""
from __future__ import annotations
import json
import structure_tools as stt

JS_URL = "https://cdn.jsdelivr.net/npm/3dmol@2.4.2/build/3Dmol-min.js"


def prepare_payload(receptor_pdb: str, pdbqt_poses_text: str, grid: dict | None = None,
                    reference_pdb_block: str | None = None, contact_cutoff: float = 4.0,
                    max_poses: int = 9, title: str = "", pdb_id: str = "") -> dict:
    poses = stt.split_pdbqt_poses(pdbqt_poses_text)[:max_poses]
    if not poses:
        raise ValueError("no docked poses found in the PDBQT text")
    info = stt.pdb_title(receptor_pdb)
    out_poses = []
    for i, p in enumerate(poses):
        res, hb = stt.ligand_contacts(receptor_pdb, p["pdb"], cutoff=contact_cutoff)
        out_poses.append({"name": f"Pose {i + 1}" + (f"  ({p['energy']:.2f} kcal/mol)" if p["energy"] is not None else ""),
                          "pdb": p["pdb"], "contacts": res, "hbonds": hb})
    return {"receptor": receptor_pdb, "poses": out_poses, "grid": grid, "reference": reference_pdb_block,
            "title": title or info["molecule"] or info["classification"] or "Receptor",
            "pdb_id": pdb_id or info["pdb_id"], "subtitle": info["title"][:140]}


_TEMPLATE = r"""
<div style="font-family:system-ui,Arial,sans-serif;font-size:13px;color:#0f172a">
<div id="hdr" style="margin-bottom:6px"></div>
<div style="display:flex;flex-wrap:wrap;gap:10px;align-items:center;margin-bottom:6px">
 <label>Pose <select id="pose"></select></label>
 <label><input type="checkbox" id="cb_grid" checked> grid box</label>
 <label><input type="checkbox" id="cb_res" checked> pocket residues</label>
 <label><input type="checkbox" id="cb_dist" checked> distances (Å)</label>
 <label id="lb_ref"><input type="checkbox" id="cb_ref" checked> crystal ligand (orange)</label>
 <label><input type="checkbox" id="cb_spin"> spin</label>
 <button id="reset">reset view</button>
</div>
<div style="position:relative;width:100%;height:520px;border:1px solid #cbd5e1;border-radius:8px">
 <div id="v" style="width:100%;height:100%;position:relative"></div>
 <div id="err" style="display:none;position:absolute;inset:0;padding:20px;background:#fff7ed;color:#9a3412">
  The 3D library could not load (blocked network?). The tables below still show every contact and distance.</div>
</div>
<div id="legend" style="margin-top:6px;color:#475569"></div>
</div>
<script src="__JS_URL__" onerror="document.getElementById('err').style.display='block'"></script>
<script>
const D = __DATA__;
(function(){
 if (typeof $3Dmol === 'undefined') { document.getElementById('err').style.display='block'; return; }
 const hdr = document.getElementById('hdr');
 hdr.innerHTML = '<b>' + D.title + '</b>' + (D.pdb_id ? ' &nbsp;<code>PDB ' + D.pdb_id + '</code>' : '') +
   (D.subtitle ? '<br><span style="color:#64748b">' + D.subtitle + '</span>' : '');
 const v = $3Dmol.createViewer(document.getElementById('v'), {backgroundColor:'white'});
 v.addModel(D.receptor, 'pdb');                      // model 0
 let refIdx = -1;
 if (D.reference) { v.addModel(D.reference, 'pdb'); refIdx = 1; } else { document.getElementById('lb_ref').style.display='none'; }
 const base = refIdx >= 0 ? 2 : 1;
 D.poses.forEach(p => v.addModel(p.pdb, 'pdb'));
 const sel = document.getElementById('pose');
 D.poses.forEach((p,i)=>{ const o=document.createElement('option'); o.value=i; o.text=p.name; sel.appendChild(o); });
 const $ = id => document.getElementById(id);
 function draw(){
  const k = parseInt(sel.value), P = D.poses[k];
  v.removeAllLabels(); v.removeAllShapes();
  v.setStyle({model:0}, {cartoon:{color:'spectrum', opacity:0.55}});
  if (refIdx >= 0) v.setStyle({model:refIdx}, $('cb_ref').checked ? {stick:{radius:0.14, color:'orange'}} : {});
  D.poses.forEach((q,i)=> v.setStyle({model:base+i}, i===k ? {stick:{radius:0.22, colorscheme:'greenCarbon'}} : {}));
  if ($('cb_res').checked) P.contacts.forEach(c=>{
    v.addStyle({model:0, chain:c.chain, resi:c.resseq}, {stick:{radius:0.12, colorscheme:'cyanCarbon'}});
    v.addLabel(c.resname + c.resseq, {position:{x:c.xyz[0],y:c.xyz[1],z:c.xyz[2]}, fontSize:11, fontColor:'#0f172a',
      backgroundColor:'white', backgroundOpacity:0.75, borderThickness:0, inFront:true});
  });
  if ($('cb_dist').checked) P.hbonds.forEach(h=>{
    v.addCylinder({start:{x:h.p_xyz[0],y:h.p_xyz[1],z:h.p_xyz[2]}, end:{x:h.l_xyz[0],y:h.l_xyz[1],z:h.l_xyz[2]},
      radius:0.07, color:'#f59e0b', dashed:true, fromCap:1, toCap:1});
    v.addLabel(h.distance_A.toFixed(2) + ' Å', {position:{x:(h.p_xyz[0]+h.l_xyz[0])/2, y:(h.p_xyz[1]+h.l_xyz[1])/2, z:(h.p_xyz[2]+h.l_xyz[2])/2},
      fontSize:12, fontColor:'#b45309', backgroundColor:'white', backgroundOpacity:0.85, borderThickness:0, inFront:true});
  });
  if (D.grid && $('cb_grid').checked) v.addBox({center:D.grid.center, dimensions:{w:D.grid.size.x,h:D.grid.size.y,d:D.grid.size.z},
      color:'#0891b2', opacity:0.0, wireframe:true});
  $('legend').innerHTML = 'Green = docked ligand &nbsp;|&nbsp; cyan sticks = pocket residues (heavy atom ≤ 4.0 Å) &nbsp;|&nbsp; ' +
    'dashed amber = N/O pairs ≤ 3.5 Å (H-bond <i>candidates</i>, distance only) &nbsp;|&nbsp; cyan box = Vina search space' +
    (refIdx >= 0 ? ' &nbsp;|&nbsp; orange = crystal pose' : '') + '<br>' + P.contacts.length + ' pocket residues, ' + P.hbonds.length + ' candidate H-bonds in ' + P.name + '.';
  v.render();
 }
 sel.onchange = draw; ['cb_grid','cb_res','cb_dist','cb_ref'].forEach(i=>$(i).onchange = draw);
 $('cb_spin').onchange = ()=> v.spin($('cb_spin').checked);
 $('reset').onclick = ()=>{ v.zoomTo({model: base + parseInt(sel.value)}); v.render(); };
 draw(); v.zoomTo({model: base}); v.zoom(0.6); v.render();
})();
</script>
"""


def build_viewer_html(receptor_pdb: str, pdbqt_poses_text: str, grid: dict | None = None,
                      reference_pdb_block: str | None = None, title: str = "", pdb_id: str = "") -> str:
    payload = prepare_payload(receptor_pdb, pdbqt_poses_text, grid, reference_pdb_block, title=title, pdb_id=pdb_id)
    data = json.dumps(payload).replace("</", "<\\/")
    return _TEMPLATE.replace("__JS_URL__", JS_URL).replace("__DATA__", data)


def show_html(html: str, height: int = 700):
    """Render HTML in Streamlit; use st.iframe when available (components.html is deprecated)."""
    import streamlit as st
    if hasattr(st, "iframe"):
        try:
            st.iframe(html, height=height)
            return
        except Exception:
            pass
    import streamlit.components.v1 as components
    components.html(html, height=height, scrolling=False)
