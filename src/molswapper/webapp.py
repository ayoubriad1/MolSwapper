"""
MolSwapper Web Interface
========================
Local web app with interactive 3D molecule visualization.

Usage:
    py -3.13 webapp.py
    Then open http://localhost:5000 in your browser.
"""

import os
import sys

from flask import Flask, render_template_string, request, jsonify

from rdkit import Chem
from rdkit.Chem import AllChem, rdDistGeom

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from molswapper.parser import parse_molecules
from molswapper.fragmenter import fragment_by_bonds, fragment_by_brics
from molswapper.recombiner import bond_cut_recombine, brics_recombine, self_recombine
from molswapper.properties import enrich_results, deduplicate
from molswapper.combinator import calculate_bond_combinations

app = Flask(__name__)

HTML_TEMPLATE = r"""
<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>MolSwapper</title>
<script src="https://3Dmol.org/build/3Dmol-min.js"></script>
<style>
  * { margin: 0; padding: 0; box-sizing: border-box; }

  body {
    font-family: 'Segoe UI', system-ui, -apple-system, sans-serif;
    background: #f0f2f5;
    color: #1a1a2e;
    min-height: 100vh;
  }

  .header {
    background: linear-gradient(135deg, #1a1a2e 0%, #16213e 100%);
    padding: 24px 32px;
    color: white;
    text-align: center;
  }

  .header h1 {
    font-size: 26px;
    font-weight: 700;
    margin-bottom: 6px;
  }

  .header .subtitle {
    font-size: 21px;
    font-weight: 600;
    opacity: 0.98;
  }

  .header .tagline {
    font-size: 13px;
    opacity: 0.65;
    margin-top: 4px;
  }

  .container {
    max-width: 1440px;
    margin: 0 auto;
    padding: 24px;
  }

  .input-panel {
    background: white;
    border: 1px solid #ddd;
    border-radius: 12px;
    padding: 24px;
    margin-bottom: 24px;
    box-shadow: 0 2px 8px rgba(0,0,0,0.06);
  }

  .input-panel h2 {
    font-size: 14px;
    font-weight: 600;
    margin-bottom: 14px;
    color: #555;
    text-transform: uppercase;
    letter-spacing: 0.5px;
  }

  .input-row {
    display: flex;
    gap: 16px;
    align-items: flex-end;
    flex-wrap: wrap;
  }

  .input-group {
    flex: 1;
    min-width: 250px;
  }

  .input-group label {
    display: block;
    font-size: 13px;
    color: #666;
    margin-bottom: 6px;
    font-weight: 500;
  }

  .input-group input, .input-group select {
    width: 100%;
    padding: 10px 14px;
    background: #f8f9fa;
    border: 1px solid #ddd;
    border-radius: 8px;
    color: #1a1a2e;
    font-size: 14px;
    font-family: 'Consolas', 'Monaco', monospace;
    outline: none;
    transition: border-color 0.2s;
  }

  .input-group select {
    font-family: 'Segoe UI', system-ui, sans-serif;
    cursor: pointer;
  }

  .input-group input:focus, .input-group select:focus {
    border-color: #4361ee;
    box-shadow: 0 0 0 3px rgba(67,97,238,0.1);
  }

  .btn {
    padding: 10px 28px;
    border: none;
    border-radius: 8px;
    font-size: 14px;
    font-weight: 600;
    cursor: pointer;
    transition: all 0.2s;
    white-space: nowrap;
  }

  .btn-primary {
    background: #4361ee;
    color: white;
  }

  .btn-primary:hover { background: #3a56d4; }

  .btn-primary:disabled {
    opacity: 0.5;
    cursor: wait;
  }

  .stats-bar {
    display: flex;
    gap: 14px;
    margin-bottom: 24px;
    flex-wrap: wrap;
  }

  .stat-card {
    background: white;
    border: 1px solid #ddd;
    border-radius: 10px;
    padding: 14px 20px;
    flex: 1;
    min-width: 140px;
    box-shadow: 0 1px 4px rgba(0,0,0,0.04);
  }

  .stat-card .label {
    font-size: 11px;
    color: #888;
    text-transform: uppercase;
    letter-spacing: 0.5px;
    margin-bottom: 4px;
  }

  .stat-card .value {
    font-size: 22px;
    font-weight: 700;
    color: #4361ee;
  }

  .section-title {
    font-size: 16px;
    font-weight: 600;
    margin: 28px 0 14px 0;
    padding-bottom: 8px;
    border-bottom: 2px solid #e0e0e0;
    display: flex;
    align-items: center;
    gap: 8px;
    color: #1a1a2e;
  }

  .section-title .badge {
    font-size: 12px;
    background: #4361ee;
    color: white;
    padding: 2px 10px;
    border-radius: 10px;
    font-weight: 600;
  }

  .parents-grid {
    display: grid;
    grid-template-columns: 1fr 1fr;
    gap: 16px;
    margin-bottom: 24px;
  }

  .parent-card {
    background: white;
    border: 1px solid #ddd;
    border-radius: 12px;
    overflow: hidden;
    box-shadow: 0 2px 8px rgba(0,0,0,0.06);
  }

  .parent-card .card-header {
    padding: 12px 16px;
    background: #f8f9fa;
    border-bottom: 1px solid #eee;
    font-weight: 600;
    font-size: 14px;
  }

  .parent-card .card-header .smi {
    font-family: 'Consolas', monospace;
    font-size: 12px;
    color: #888;
    font-weight: 400;
    margin-top: 4px;
    word-break: break-all;
  }

  .mol-viewer {
    width: 100%;
    height: 300px;
    position: relative;
    background: black;
    border-bottom: 1px solid #eee;
  }

  .viewer-loading {
    position: absolute;
    top: 50%;
    left: 50%;
    transform: translate(-50%, -50%);
    text-align: center;
    color: #aaa;
    font-size: 13px;
    z-index: 1;
    pointer-events: none;
  }

  .viewer-loading .mini-spinner {
    width: 24px;
    height: 24px;
    border: 3px solid #e0e0e0;
    border-top-color: #4361ee;
    border-radius: 50%;
    animation: spin 0.7s linear infinite;
    margin: 0 auto 8px;
  }

  @keyframes spin { to { transform: rotate(360deg); } }

  .fragments-list {
    padding: 12px 16px;
  }

  .fragments-list .frag-title {
    font-size: 12px;
    color: #888;
    text-transform: uppercase;
    letter-spacing: 0.5px;
    margin-bottom: 6px;
  }

  .frag-item {
    display: flex;
    align-items: center;
    gap: 8px;
    padding: 3px 0;
    font-size: 13px;
  }

  .frag-dot {
    width: 8px;
    height: 8px;
    border-radius: 50%;
    flex-shrink: 0;
  }

  .frag-dot.a { background: #4361ee; }
  .frag-dot.b { background: #22c55e; }

  .frag-smi {
    font-family: 'Consolas', monospace;
    font-size: 12px;
    color: #888;
  }

  .products-grid {
    display: grid;
    grid-template-columns: repeat(auto-fill, minmax(350px, 1fr));
    gap: 16px;
  }

  .product-card {
    background: white;
    border: 1px solid #ddd;
    border-radius: 12px;
    overflow: hidden;
    box-shadow: 0 1px 6px rgba(0,0,0,0.05);
    transition: box-shadow 0.2s, border-color 0.2s;
  }

  .product-card:hover {
    border-color: #4361ee;
    box-shadow: 0 4px 16px rgba(67,97,238,0.12);
  }

  .product-card .card-header {
    padding: 10px 16px;
    background: #f8f9fa;
    border-bottom: 1px solid #eee;
  }

  .product-card .card-header .pid {
    font-size: 12px;
    font-weight: 600;
    color: #333;
  }

  .product-card .card-header .origin {
    font-size: 12px;
    color: #666;
    margin-top: 3px;
    display: flex;
    align-items: center;
    gap: 6px;
  }

  .origin-tag {
    display: inline-block;
    padding: 2px 8px;
    border-radius: 4px;
    font-size: 11px;
    font-weight: 600;
  }

  .origin-tag.a {
    background: rgba(67, 97, 238, 0.12);
    color: #4361ee;
  }

  .origin-tag.b {
    background: rgba(34, 197, 94, 0.12);
    color: #16a34a;
  }

  .product-viewer {
    height: 250px;
    background: black;
    position: relative;
  }

  .product-props {
    padding: 8px 16px;
    border-top: 1px solid #eee;
    display: grid;
    grid-template-columns: repeat(4, 1fr);
    gap: 4px;
  }

  .prop {
    text-align: center;
    padding: 4px;
  }

  .prop .prop-label {
    font-size: 10px;
    color: #999;
    text-transform: uppercase;
  }

  .prop .prop-value {
    font-size: 14px;
    font-weight: 600;
    color: #333;
  }

  .product-smiles {
    padding: 8px 16px;
    border-top: 1px solid #eee;
    font-family: 'Consolas', monospace;
    font-size: 11px;
    color: #999;
    word-break: break-all;
    cursor: pointer;
    transition: color 0.15s;
  }

  .product-smiles:hover { color: #4361ee; }

  .swap-section { margin-bottom: 28px; }

  .swap-label {
    font-size: 14px;
    font-weight: 600;
    margin-bottom: 12px;
    padding: 8px 16px;
    background: white;
    border: 1px solid #ddd;
    border-radius: 8px;
    display: inline-flex;
    align-items: center;
    gap: 8px;
    color: #333;
  }

  .loading-state {
    text-align: center;
    padding: 60px;
    color: #888;
  }

  .loading-state .spinner {
    width: 36px;
    height: 36px;
    border: 4px solid #e0e0e0;
    border-top-color: #4361ee;
    border-radius: 50%;
    animation: spin 0.7s linear infinite;
    margin: 0 auto 16px;
  }

  #results { display: none; }

  .empty-state {
    text-align: center;
    padding: 80px 24px;
    color: #888;
  }

  .empty-state h3 {
    font-size: 18px;
    margin-bottom: 8px;
    color: #333;
  }

  .example-btn {
    background: white;
    border: 1px solid #ddd;
    color: #4361ee;
    padding: 8px 18px;
    border-radius: 6px;
    cursor: pointer;
    font-size: 13px;
    font-weight: 500;
    margin-top: 16px;
    transition: all 0.2s;
  }

  .example-btn:hover {
    background: #f0f4ff;
    border-color: #4361ee;
  }

  .progress-bar-wrap {
    width: 100%;
    background: #e0e0e0;
    border-radius: 6px;
    height: 6px;
    margin-top: 12px;
    overflow: hidden;
  }

  .progress-bar-fill {
    height: 100%;
    background: #4361ee;
    border-radius: 6px;
    transition: width 0.3s;
    width: 0%;
  }

  @media (max-width: 768px) {
    .parents-grid { grid-template-columns: 1fr; }
    .products-grid { grid-template-columns: 1fr; }
    .input-row { flex-direction: column; }
    .product-props { grid-template-columns: repeat(2, 1fr); }
  }
</style>
</head>
<body>

<div class="header">
  <h1>MolSwapper</h1>
  <div class="subtitle">Made by Ayoub Riad &mdash; Researcher in Bioinformatics</div>
  <div class="tagline">Molecular Fragment Swapping with 3D Visualization</div>
</div>

<div class="container">

  <div class="input-panel">
    <h2>Input Molecules</h2>
    <div class="input-row">
      <div class="input-group">
        <label>Molecule A (SMILES)</label>
        <input type="text" id="smilesA" placeholder="e.g. c1ccccc1CCO">
      </div>
      <div class="input-group">
        <label>Molecule B (SMILES)</label>
        <input type="text" id="smilesB" placeholder="e.g. c1ccccc1CCNC">
      </div>
      <div class="input-group" style="flex:0.5; min-width: 140px;">
        <label>Method</label>
        <select id="method">
          <option value="bonds">Bond-cut</option>
          <option value="brics">BRICS</option>
          <option value="both">Both</option>
        </select>
      </div>
      <button class="btn btn-primary" id="runBtn" onclick="runSwap()">
        Fragment &amp; Swap
      </button>
    </div>
  </div>

  <div id="emptyState" class="empty-state">
    <h3>Enter two molecules to begin</h3>
    <p>Paste SMILES strings above, then click Fragment &amp; Swap to generate chimeric molecules.</p>
    <button class="example-btn" onclick="loadExample()">Load example molecules</button>
  </div>

  <div id="loadingState" class="loading-state" style="display:none;">
    <div class="spinner"></div>
    <p>Fragmenting and recombining molecules...</p>
    <p style="font-size:12px; margin-top:4px;">Generating 3D conformers on the server</p>
  </div>

  <div id="results">
    <div class="stats-bar" id="statsBar"></div>

    <div class="section-title">Input Molecules</div>
    <div class="parents-grid" id="parentsGrid"></div>

    <div class="section-title">
      Generated Chimeras
      <span class="badge" id="productCount">0</span>
    </div>
    <div id="progressWrap" style="display:none;">
      <div style="font-size:13px; color:#888; margin-bottom:4px;">Rendering 3D viewers...</div>
      <div class="progress-bar-wrap"><div class="progress-bar-fill" id="progressBar"></div></div>
    </div>
    <div id="productsContainer" style="margin-top:16px;"></div>
  </div>
</div>

<script>
function loadExample() {
  document.getElementById('smilesA').value = 'c1ccccc1CCO';
  document.getElementById('smilesB').value = 'c1ccccc1CCNC';
}

function copySmiles(el) {
  navigator.clipboard.writeText(el.textContent.trim());
  const orig = el.style.color;
  el.style.color = '#22c55e';
  setTimeout(() => { el.style.color = orig; }, 600);
}

function createViewer(container, molBlock) {
  if (!molBlock) {
    container.innerHTML = '<div style="display:flex;align-items:center;justify-content:center;height:100%;color:#ccc;font-size:13px;">No 3D data</div>';
    return;
  }

  try {
    container.id = 'viewer_' + Math.random().toString(36).substr(2, 9);
    var viewer = $3Dmol.createViewer(container.id, {
      backgroundColor: 'black',
      antialias: true,
      disableFog: true
    });

    viewer.addModel(molBlock, 'sdf');
    viewer.setStyle({}, {
      sphere: { scale: 0.3, colorscheme: 'greenCarbon' },
      stick:  { radius: 0.15, colorscheme: 'greenCarbon' }
    });
    viewer.zoomTo();
    viewer.zoom(0.8);
    viewer.render();
  } catch(e) {
    console.error('3Dmol viewer error:', e);
    container.innerHTML = '<div style="display:flex;align-items:center;justify-content:center;height:100%;color:#e44;font-size:13px;">Viewer error</div>';
  }
}

async function runSwap() {
  var smiA = document.getElementById('smilesA').value.trim();
  var smiB = document.getElementById('smilesB').value.trim();
  var method = document.getElementById('method').value;

  if (!smiA || !smiB) {
    alert('Please enter SMILES for both molecules.');
    return;
  }

  document.getElementById('emptyState').style.display = 'none';
  document.getElementById('results').style.display = 'none';
  document.getElementById('loadingState').style.display = 'block';
  document.getElementById('runBtn').disabled = true;

  try {
    var resp = await fetch('/api/swap', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ smiles_a: smiA, smiles_b: smiB, method: method })
    });

    var data = await resp.json();
    if (data.error) {
      alert('Error: ' + data.error);
      document.getElementById('loadingState').style.display = 'none';
      document.getElementById('emptyState').style.display = 'block';
      return;
    }

    renderResults(data);
  } catch (e) {
    alert('Request failed: ' + e.message);
    document.getElementById('loadingState').style.display = 'none';
    document.getElementById('emptyState').style.display = 'block';
  } finally {
    document.getElementById('runBtn').disabled = false;
  }
}

function renderResults(data) {
  document.getElementById('loadingState').style.display = 'none';
  document.getElementById('results').style.display = 'block';

  // Stats bar
  var stats = data.stats;
  document.getElementById('statsBar').innerHTML =
    '<div class="stat-card"><div class="label">Theoretical Max</div><div class="value">' + stats.theoretical_max + '</div></div>' +
    '<div class="stat-card"><div class="label">Valid Products</div><div class="value">' + stats.total_products + '</div></div>' +
    '<div class="stat-card"><div class="label">Unique</div><div class="value">' + stats.unique_products + '</div></div>' +
    '<div class="stat-card"><div class="label">MW Range</div><div class="value" style="font-size:15px">' + stats.mw_min + ' - ' + stats.mw_max + '</div></div>';

  // Build parent cards using DOM (not innerHTML) so viewers stay alive
  var pg = document.getElementById('parentsGrid');
  pg.innerHTML = '';

  data.parents.forEach(function(p, i) {
    var tag = (i === 0) ? 'A' : 'B';
    var cls = (i === 0) ? 'a' : 'b';

    var card = document.createElement('div');
    card.className = 'parent-card';

    var header = document.createElement('div');
    header.className = 'card-header';
    header.innerHTML =
      '<span class="origin-tag ' + cls + '">Molecule ' + tag + '</span> ' + p.name +
      '<div class="smi">' + p.smiles + '</div>';
    card.appendChild(header);

    var viewerDiv = document.createElement('div');
    viewerDiv.className = 'mol-viewer';
    card.appendChild(viewerDiv);

    var fragDiv = document.createElement('div');
    fragDiv.className = 'fragments-list';
    var fragHtml = '<div class="frag-title">Fragments (' + p.fragments.length + ' cuts)</div>';
    p.fragments.forEach(function(f, fi) {
      fragHtml += '<div class="frag-item"><span class="frag-dot ' + cls + '"></span>Cut ' + (fi+1) +
        ': <span class="frag-smi">' + f[0] + '</span> + <span class="frag-smi">' + f[1] + '</span></div>';
    });
    fragDiv.innerHTML = fragHtml;
    card.appendChild(fragDiv);

    pg.appendChild(card);

    // Init viewer after it is in the DOM
    setTimeout(function() { createViewer(viewerDiv, p.mol_block); }, 100);
  });

  // Products - build progressively, one at a time
  document.getElementById('productCount').textContent = data.products.length;
  var container = document.getElementById('productsContainer');
  container.innerHTML = '';

  // Group by cut_pair
  var cutGroups = {};
  var groupOrder = [];
  data.products.forEach(function(p) {
    var key = p.cut_pair || 'brics';
    if (!cutGroups[key]) {
      cutGroups[key] = [];
      groupOrder.push(key);
    }
    cutGroups[key].push(p);
  });

  // Build all cards first (without viewers), then init viewers one-by-one
  var allViewerJobs = [];

  groupOrder.forEach(function(cutPair) {
    var products = cutGroups[cutPair];

    var section = document.createElement('div');
    section.className = 'swap-section';

    var label = document.createElement('div');
    label.className = 'swap-label';
    label.textContent = (cutPair === 'brics') ? 'BRICS Recombination' : cutPair;
    section.appendChild(label);

    var grid = document.createElement('div');
    grid.className = 'products-grid';

    products.forEach(function(p) {
      var card = document.createElement('div');
      card.className = 'product-card';

      var hdr = document.createElement('div');
      hdr.className = 'card-header';
      hdr.innerHTML =
        '<div class="pid">' + p.product_id + '</div>' +
        '<div class="origin">' +
          '<span class="origin-tag a">' + p.frag_from_a + '</span>' +
          '<span style="color:#aaa;">+</span>' +
          '<span class="origin-tag b">' + p.frag_from_b + '</span>' +
        '</div>';
      card.appendChild(hdr);

      var viewerDiv = document.createElement('div');
      viewerDiv.className = 'product-viewer';
      viewerDiv.innerHTML = '<div class="viewer-loading"><div class="mini-spinner"></div>Generating 3D...</div>';
      card.appendChild(viewerDiv);

      var props = document.createElement('div');
      props.className = 'product-props';
      props.innerHTML =
        '<div class="prop"><div class="prop-label">MW</div><div class="prop-value">' + p.mol_weight + '</div></div>' +
        '<div class="prop"><div class="prop-label">LogP</div><div class="prop-value">' + p.logp + '</div></div>' +
        '<div class="prop"><div class="prop-label">HBD</div><div class="prop-value">' + p.hbd + '</div></div>' +
        '<div class="prop"><div class="prop-label">HBA</div><div class="prop-value">' + p.hba + '</div></div>';
      card.appendChild(props);

      var smiDiv = document.createElement('div');
      smiDiv.className = 'product-smiles';
      smiDiv.title = 'Click to copy';
      smiDiv.textContent = p.product_smiles;
      smiDiv.onclick = function() { copySmiles(this); };
      card.appendChild(smiDiv);

      grid.appendChild(card);

      allViewerJobs.push({ el: viewerDiv, block: p.mol_block });
    });

    section.appendChild(grid);
    container.appendChild(section);
  });

  // Progressive viewer initialization: render them one by one
  var progressWrap = document.getElementById('progressWrap');
  var progressBar = document.getElementById('progressBar');
  var total = allViewerJobs.length;

  if (total > 0) {
    progressWrap.style.display = 'block';
    var idx = 0;

    function renderNext() {
      if (idx >= total) {
        progressWrap.style.display = 'none';
        return;
      }
      try {
        var job = allViewerJobs[idx];
        job.el.innerHTML = '';
        createViewer(job.el, job.block);
      } catch(e) {
        console.error('Render error at index ' + idx, e);
      }
      idx++;
      progressBar.style.width = Math.round((idx / total) * 100) + '%';
      setTimeout(renderNext, 150);
    }

    setTimeout(renderNext, 300);
  }
}

document.addEventListener('keydown', function(e) {
  if (e.key === 'Enter' && !document.getElementById('runBtn').disabled) {
    runSwap();
  }
});
</script>

</body>
</html>
"""


def _mol_to_block(smiles: str) -> str:
    """Generate a 3D mol block from SMILES for the 3Dmol.js viewer."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return ""
    mol = Chem.AddHs(mol)
    params = rdDistGeom.ETKDGv3()
    params.randomSeed = 42
    params.useSmallRingTorsions = True
    params.numThreads = 0
    cids = rdDistGeom.EmbedMultipleConfs(mol, numConfs=1, params=params)
    if len(cids) == 0:
        params.useRandomCoords = True
        cids = rdDistGeom.EmbedMultipleConfs(mol, numConfs=1, params=params)
        if len(cids) == 0:
            return ""
    try:
        props = AllChem.MMFFGetMoleculeProperties(mol)
        if props:
            ff = AllChem.MMFFGetMoleculeForceField(mol, props, confId=cids[0])
            if ff:
                ff.Minimize(maxIts=500)
    except Exception:
        pass
    return Chem.MolToMolBlock(mol, confId=cids[0])


@app.route("/")
def index():
    return render_template_string(HTML_TEMPLATE)


@app.route("/api/swap", methods=["POST"])
def api_swap():
    data = request.get_json()
    smi_a = data.get("smiles_a", "").strip()
    smi_b = data.get("smiles_b", "").strip()
    method = data.get("method", "bonds")

    mol_a = Chem.MolFromSmiles(smi_a)
    mol_b = Chem.MolFromSmiles(smi_b)
    if mol_a is None:
        return jsonify({"error": f"Invalid SMILES for Molecule A: {smi_a}"})
    if mol_b is None:
        return jsonify({"error": f"Invalid SMILES for Molecule B: {smi_b}"})

    can_a = Chem.MolToSmiles(mol_a)
    can_b = Chem.MolToSmiles(mol_b)

    mol_data = {"mol_A": mol_a, "mol_B": mol_b}

    frags_a = fragment_by_bonds(mol_a)
    frags_b = fragment_by_bonds(mol_b)

    combo = calculate_bond_combinations(mol_data)

    all_results = []
    if method in ("bonds", "both"):
        all_results.extend(bond_cut_recombine(mol_data))
    if method in ("brics", "both"):
        all_results.extend(brics_recombine(mol_data, max_per_pair=50))

    all_results = enrich_results(all_results)
    unique_results = deduplicate(all_results)

    block_a = _mol_to_block(can_a)
    block_b = _mol_to_block(can_b)

    products_out = []
    for r in unique_results:
        smi = r["product_smiles"]
        mb = _mol_to_block(smi)

        swap = r.get("swap_type", "")

        if swap in ("A1-B1", "A1-B2", "A2-B1", "A2-B2"):
            side_a = "Left" if swap[1] == "1" else "Right"
            side_b = "Left" if swap[4] == "1" else "Right"
            frag_from_a = f"{side_a} of A"
            frag_from_b = f"{side_b} of B"
        else:
            frag_from_a = "BRICS frag A"
            frag_from_b = "BRICS frag B"

        pid = r.get("product_id", "")
        if "cut" in pid and "_mol_" in pid:
            parts = pid.split("_")
            cut_indices = [i for i, p in enumerate(parts) if p.startswith("cut")]
            if len(cut_indices) >= 2:
                cut_pair = f"A {parts[cut_indices[0]]} × B {parts[cut_indices[1]]}"
            else:
                cut_pair = r.get("frag_a", "unknown")
        else:
            cut_pair = "brics"

        products_out.append({
            "product_id": pid,
            "product_smiles": smi,
            "mol_block": mb,
            "swap_type": swap,
            "frag_from_a": frag_from_a,
            "frag_from_b": frag_from_b,
            "frag_a": r.get("frag_a", ""),
            "frag_b": r.get("frag_b", ""),
            "cut_pair": cut_pair,
            "mol_weight": r.get("mol_weight", 0),
            "logp": r.get("logp", 0),
            "hbd": r.get("hbd", 0),
            "hba": r.get("hba", 0),
            "tpsa": r.get("tpsa", 0),
            "num_atoms": r.get("num_atoms", 0),
        })

    weights = [p["mol_weight"] for p in products_out] if products_out else [0]
    stats = {
        "theoretical_max": combo["total_theoretical_max"],
        "total_products": len(all_results),
        "unique_products": len(unique_results),
        "mw_min": min(weights),
        "mw_max": max(weights),
    }

    return jsonify({
        "parents": [
            {
                "name": "mol_A",
                "smiles": can_a,
                "mol_block": block_a,
                "fragments": [(f[0], f[1]) for f in frags_a],
            },
            {
                "name": "mol_B",
                "smiles": can_b,
                "mol_block": block_b,
                "fragments": [(f[0], f[1]) for f in frags_b],
            },
        ],
        "products": products_out,
        "stats": stats,
    })


def main():
    print("\n  MolSwapper Web Interface")
    print("  Developed by Ayoub Riad, a Researcher in Bioinformatics")
    print("  Open http://localhost:5000 in your browser\n")
    app.run(host="127.0.0.1", port=5000, debug=False)


if __name__ == "__main__":
    main()
