"""Self-contained multi-view frontend used by the 3C interface RIN exporter."""

from __future__ import annotations

from pathlib import Path


CYTOSCAPE_VERSION = "3.34.4"
CYTOSCAPE_ASSET = Path(__file__).with_name("pymacs_visualization_assets") / f"cytoscape-{CYTOSCAPE_VERSION}.min.js"
THREEDMOL_ASSET = Path(__file__).with_name("pymacs_visualization_assets") / "3Dmol-min.js"


def _build_interface_explorer_legacy(payload: str) -> str:
    """Return the offline-capable explorer.  ``payload`` is trusted generated JSON."""
    cytoscape_source = CYTOSCAPE_ASSET.read_text(encoding="utf-8")
    return r'''<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>PyMACS Interface Explorer</title>
<style>:root{--bg:#05070d;--p:#0e1625;--l:#2a405e;--t:#ecf5ff;--m:#9bb0c9;--c:#4bd7ff;--pep:#ffcf5b}*{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--t);font:14px Inter,system-ui,sans-serif}main{max-width:1500px;margin:auto;padding:24px}.ey{color:var(--c);font-weight:900;letter-spacing:.13em;font-size:11px}h1{font-size:clamp(27px,4vw,52px);margin:5px 0}p{color:var(--m);line-height:1.5}.bar{display:flex;gap:8px;flex-wrap:wrap;margin:14px 0}button,select,input{border:1px solid var(--l);background:#101d30;color:var(--t);border-radius:8px;padding:8px;font:inherit}button{font-weight:800;cursor:pointer}button.on{border-color:var(--c);background:#12425b}.panel{background:var(--p);border:1px solid var(--l);border-radius:15px;padding:12px}.view{min-height:620px}.hidden{display:none!important}.split{display:grid;grid-template-columns:minmax(0,1fr) 310px;gap:12px}.detail{padding:12px;border-left:1px solid var(--l);color:var(--m)}svg{width:100%;height:610px}.node{cursor:pointer}.node circle{stroke:#07101e;stroke-width:2}.edge{fill:none;cursor:pointer}.seq{display:flex;flex-wrap:wrap;gap:5px}.res{padding:6px;min-width:35px;text-align:center;background:#17243a;border-radius:6px;cursor:pointer}.res.pep{border:1px solid var(--pep)}.res.sel{outline:3px solid var(--c)}table{border-collapse:collapse}td,th{padding:5px;border:1px solid #29415b;font-size:11px}td{cursor:pointer;min-width:24px;height:24px}canvas{width:100%;height:auto;background:#07111e}#space{touch-action:none}#cy{height:610px;min-width:0;background:#07111e;border-radius:10px}#structure{min-height:0}.structure-grid{display:grid;grid-template-columns:minmax(0,1fr) minmax(245px,30%);gap:12px;align-items:stretch}.structure-stage{position:relative;min-width:0;height:clamp(380px,52vh,620px);overflow:hidden;border:1px solid var(--l);border-radius:10px;background:#05070d}.structure-stage #mol,.structure-stage #mol canvas{width:100%!important;height:100%!important;display:block}.structure-detail{max-height:clamp(380px,52vh,620px);overflow:auto;padding:12px;border:1px solid var(--l);border-radius:10px;background:#0b1423;color:var(--m)}.structure-detail h3{margin:0 0 8px;color:var(--t)}.structure-detail h4{margin:14px 0 6px;color:var(--t);font-size:12px}.structure-detail .seq{gap:6px}.structure-detail .res{font-size:12px;padding:5px;min-width:32px}#rin{display:none}@media(max-width:900px){.split,.structure-grid{grid-template-columns:1fr}.detail{border-left:0;border-top:1px solid var(--l)}.structure-detail{max-height:300px}.structure-stage{height:440px}}</style>
<main><div class="ey">PYMACS · PEPTIDE–PROTEIN INTERACTION EXPLORER</div><h1>Persistent interface explorer</h1><p id="note">Geometric contact occupancy summarizes an MD ensemble. Validated interaction classes are retained separately; they are not inferred from the displayed structure.</p><div class="bar" id="tabs"><button data-v="network" class="on">Network</button><button data-v="structure">3D Structure</button><button data-v="spatial">3D Nodes</button><button data-v="map">Contact Map</button><button data-v="sequence">Peptide Sequence</button><button id="export">Export SVG</button></div><div class="bar"><label>Interface <select id="pair"></select></label><label>Minimum contact occupancy <input id="occ" type="range" min="0" max="100" step="1"><b id="occv"></b></label><label>Time window <select id="window"><option value="all">Entire trajectory</option><option value="first">First 25%</option><option value="middle">Middle 50%</option><option value="final">Final 25%</option><option value="custom">Custom ns</option></select></label><span id="custom" class="hidden"><input id="start" type="number" step="any" aria-label="Start ns">–<input id="end" type="number" step="any" aria-label="End ns"></span><label>Type <select id="type"><option value="all">All contact evidence</option><option>Hydrogen bond</option><option>Ionic</option><option>Hydrophobic</option><option>Other contact</option></select></label></div><section class="panel"><div id="network" class="view split"><div><div class="bar" id="layouts"><button data-layout="cose" class="on">Force-directed</button><button data-layout="bipartite">Bipartite</button><button data-layout="concentric">Concentric</button><button data-layout="circle">Circular</button><button data-layout="chain">Chain-separated</button></div><div id="cy"></div><svg id="rin" viewBox="0 0 1200 610"></svg></div><aside class="detail" id="detail">Select a residue or edge to inspect it.</aside></div><div id="structure" class="view hidden"><div class="structure-grid"><div class="structure-stage"><div id="mol" aria-label="3D protein-peptide structure"></div></div><aside class="structure-detail" id="structure-detail">Select a residue to inspect it.</aside></div></div><div id="spatial" class="view hidden"><div class="bar"><button id="zoom-in">Zoom in</button><button id="zoom-out">Zoom out</button><button id="zoom-reset">Fit contacts</button><span>Scroll to zoom · drag to rotate</span></div><canvas id="space" width="1200" height="610"></canvas></div><div id="map" class="view hidden"></div><div id="sequence" class="view hidden"></div></section></main>
<script>''' + cytoscape_source + r'''</script><script>const DATA=''' + payload + r''';
(()=>{let state={pair:DATA.pairs[0]?.id,occ:DATA.default_threshold||.5,kind:'all',selected:null,view:'network',a:.35},viewer;const $=s=>document.querySelector(s),esc=s=>String(s).replace(/[&<>"']/g,x=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[x])),times=DATA.presence?.times_ns||[],rows=DATA.presence?.rows||{},sep='\u001f';DATA.edges.forEach(e=>e.aggregate=e.fraction);const pair=()=>DATA.pairs.find(p=>p.id===state.pair);const key=e=>e.a+sep+e.b;function range(){let a=times[0],b=times.at(-1),w=$('#window').value;if(w==='first')b=a+(b-a)/4;if(w==='middle'){let d=b-a;a+=d/4;b-=d/4}if(w==='final')a=b-(b-a)/4;if(w==='custom'){a=+$('#start').value;b=+$('#end').value}return[a,b]}function recalc(){if(!times.length||!Object.keys(rows).length){$('#note').textContent='Time-window filtering unavailable: this historical run has no saved per-frame geometric-contact matrix. Aggregate occupancy is shown.';return}let[a,b]=range();DATA.edges.forEach(e=>{let bin=rows[key(e)]||rows[e.b+sep+e.a];if(!bin)return;bin=Uint8Array.from(atob(bin),c=>c.charCodeAt(0));let n=0,c=0;times.forEach((t,i)=>{if(t>=a&&t<=b){n++;c+=bin[i]||0}});e.fraction=n?c/n:0});$('#note').textContent=`Geometric contact occupancy: ${a.toFixed(3)}–${b.toFixed(3)} ns. Chemistry labels are validated-event evidence, never structure-derived.`}function edges(){let p=pair();return DATA.edges.filter(e=>p&&e.a.startsWith(p.left+':')&&e.b.startsWith(p.right+':')&&e.fraction>=state.occ&&(state.kind==='all'||e.type===state.kind))}function select(id,e){state.selected=id;let ps=edges().filter(x=>x.a===id||x.b===id);$('#detail').innerHTML='<b>'+esc(id)+'</b><br>'+ (ps.length?ps.map(x=>esc(x.a===id?x.b:x.a)+' · '+(x.fraction*100).toFixed(1)+'% · '+esc(x.type)).join('<br>'):'No partners under active filters.');window.parent?.postMessage({type:'pymacs-rin-selection',residue:id,edge:e&&[e.a,e.b]},'*');draw()}function draw(){let p=pair(),es=edges(),svg=$('#rin');if(!p)return;let A=[...new Set(es.map(e=>e.a))],B=[...new Set(es.map(e=>e.b))],pos={},at=(x,i,n)=>[x,80+i*470/Math.max(n-1,1)];A.forEach((x,i)=>pos[x]=at(180,i,A.length));B.forEach((x,i)=>pos[x]=at(1000,i,B.length));svg.innerHTML=es.map((e,i)=>{let a=pos[e.a],b=pos[e.b],col=/bond/i.test(e.type)?'#ff5578':/ionic/i.test(e.type)?'#c084fc':/hydrophobic/i.test(e.type)?'#7dff91':'#94a3b8';return `<path class="edge" data-e="${i}" d="M${a} Q600 ${(a[1]+b[1])/2} ${b}" stroke="${col}" stroke-width="${2+7*e.fraction}" opacity="${.3+.7*e.fraction}"/>`}).join('')+[...A,...B].map(id=>{let q=pos[id],sel=state.selected===id;return `<g class="node" data-r="${esc(id)}"><circle cx="${q[0]}" cy="${q[1]}" r="${sel?32:25}" fill="${A.includes(id)?'#ffcf5b':'#4bd7ff'}"/><text x="${q[0]}" y="${q[1]+4}" text-anchor="middle" fill="#07101e" font-weight="900">${esc(id.split(':').at(-1))}</text></g>`}).join('');svg.querySelectorAll('[data-r]').forEach(n=>n.onclick=()=>select(n.dataset.r));svg.querySelectorAll('[data-e]').forEach(n=>{let e=es[n.dataset.e];n.onclick=()=>select(e.a,e)})}function sequence(){let p=pair();$('#sequence').innerHTML=[p.left,p.right].map(c=>{let ch=DATA.chains[c];return `<h3>${esc(c)} · ${esc(ch.type)}</h3><div class="seq">${ch.residues.map(r=>`<span class="res ${ch.type==='peptide'?'pep':''} ${state.selected===r.label?'sel':''}" data-r="${esc(r.label)}" title="${esc(r.stable_id||r.label)}">${esc(r.one)}${r.number}</span>`).join('')}</div>`}).join('');$('#sequence').querySelectorAll('[data-r]').forEach(n=>n.onclick=()=>select(n.dataset.r))}function contact(){let p=pair(),L=DATA.chains[p.left].residues,R=DATA.chains[p.right].residues,by=new Map(edges().flatMap(e=>[[key(e),e],[e.b+sep+e.a,e]]));$('#map').innerHTML='<p>Rows preserve peptide N→C order when a peptide is the left selected chain. Cell colour is geometric-contact occupancy.</p><table><tr><th></th>'+R.map(r=>'<th>'+esc(r.one+r.number)+'</th>').join('')+'</tr>'+L.map(a=>'<tr><th>'+esc(a.one+a.number)+'</th>'+R.map(b=>{let e=by.get(a.label+sep+b.label),v=e?.fraction||0;return `<td data-a="${esc(a.label)}" title="${esc(a.label)} ↔ ${esc(b.label)} · ${(100*v).toFixed(2)}%" style="background:rgba(75,215,255,${.05+.9*v})"></td>`}).join('')+'</tr>').join('')+'</table>';$('#map').querySelectorAll('td').forEach(n=>n.onclick=()=>select(n.dataset.a))}function spatial(){let p=pair(),all=[...DATA.chains[p.left].residues,...DATA.chains[p.right].residues],by=new Map(all.map(r=>[r.label,r])),c=$('#space').getContext('2d'),co=r=>r.position||[0,0,0],pr=r=>{let[x,y,z]=co(r),q=Math.cos(state.a),s=Math.sin(state.a);return[600+(x*q-z*s)*7,305-y*7+(x*s+z*q)]};c.clearRect(0,0,1200,610);edges().forEach(e=>{let a=by.get(e.a),b=by.get(e.b);if(!a||!b)return;let A=pr(a),B=pr(b);c.strokeStyle='#7bdfff';c.globalAlpha=.25+.75*e.fraction;c.lineWidth=1+5*e.fraction;c.beginPath();c.moveTo(...A);c.lineTo(...B);c.stroke()});c.globalAlpha=1;all.forEach(r=>{let[x,y]=pr(r),isPep=DATA.chains[p.left].type==='peptide'&&DATA.chains[p.left].residues.includes(r)||DATA.chains[p.right].type==='peptide'&&DATA.chains[p.right].residues.includes(r);c.fillStyle=isPep?'#ffcf5b':'#4bd7ff';c.beginPath();c.arc(x,y,state.selected===r.label?9:5,0,7);c.fill();c.fillStyle='#fff';c.fillText(r.one+r.number,x+7,y-7)})}function structure(){let host=$('#structure');if(viewer)return;if(!DATA.structure_pdb_b64){host.textContent='No coordinate PDB was retained for this run.';return}let load=()=>{if(!window.$3Dmol){host.textContent='3Dmol.js is not locally available in this standalone report.';return}viewer=$3Dmol.createViewer(host,{backgroundColor:'#05070d'});viewer.addModel(atob(DATA.structure_pdb_b64),'pdb');viewer.setStyle({},{cartoon:{color:'#53718c'}});viewer.zoomTo();viewer.render()};let s=document.createElement('script');s.src='/static/3Dmol-min.js';s.onload=load;s.onerror=load;document.head.appendChild(s)}function render(){draw();sequence();contact();spatial();if(state.view==='structure')structure()}$('#pair').innerHTML=DATA.pairs.map(p=>`<option value="${esc(p.id)}">${esc(p.left)} ↔ ${esc(p.right)}</option>`).join('');$('#occ').value=state.occ*100;$('#occv').textContent=Math.round(state.occ*100)+'%';if(times.length){$('#start').value=times[0];$('#end').value=times.at(-1)}$('#pair').onchange=e=>{state.pair=e.target.value;render()};$('#occ').oninput=e=>{state.occ=e.target.value/100;$('#occv').textContent=e.target.value+'%';render()};$('#type').onchange=e=>{state.kind=e.target.value;render()};$('#window').onchange=e=>{$('#custom').classList.toggle('hidden',e.target.value!=='custom');recalc();render()};$('#start').onchange=$('#end').onchange=()=>{recalc();render()};$('#tabs').querySelectorAll('button[data-v]').forEach(b=>b.onclick=()=>{state.view=b.dataset.v;$('#tabs').querySelectorAll('button[data-v]').forEach(x=>x.classList.toggle('on',x===b));['network','structure','spatial','map','sequence'].forEach(v=>$('#'+v).classList.toggle('hidden',v!==state.view));render()});$('#space').onmousedown=e=>{let x=e.clientX,a=state.a;let m=q=>{state.a=a+(q.clientX-x)/150;spatial()},u=()=>{removeEventListener('mousemove',m);removeEventListener('mouseup',u)};addEventListener('mousemove',m);addEventListener('mouseup',u)};$('#export').onclick=()=>{let a=document.createElement('a'),b=new Blob([$('#rin').outerHTML],{type:'image/svg+xml'});a.href=URL.createObjectURL(b);a.download='Interface_RIN.svg';a.click()};recalc();render()})();</script></html>'''
CYTOSCAPE_BOOTSTRAP = r'''<script>
/* Cytoscape.js 3.34.4 primary renderer; the generated SVG stays as fallback. */
(()=>{if(!window.cytoscape||document.querySelector('#cy').offsetWidth===0)return;let cy,layout='bipartite';const $=s=>document.querySelector(s),allEdges=Array.isArray(DATA.edges)?DATA.edges:Object.values(DATA.edges||{});const active=()=>{let p=DATA.pairs.find(x=>x.id===$('#pair').value),o=+$('#occ').value/100,k=$('#type').value;return{p,es:allEdges.filter(e=>p&&((e.a.startsWith(p.left+':')&&e.b.startsWith(p.right+':'))||(e.a.startsWith(p.right+':')&&e.b.startsWith(p.left+':')))&&e.fraction>=o&&(k==='all'||e.type===k))}};const cls=t=>/bond/i.test(t)?'bond':/ionic/i.test(t)?'ionic':/hydrophobic/i.test(t)?'hydrophobic':'other';function run(){let{p,es}=active();if(!p)return;$('#detail').textContent=es.length?'Select a residue or edge to inspect it.':'No contacts meet the active filter. Lower the occupancy threshold to inspect the interface.';let l=DATA.chains[p.left].residues,r=DATA.chains[p.right].residues,visible=new Set(es.flatMap(e=>[e.a,e.b])),nodes=[...l,...r].filter(x=>visible.has(x.label)).map(x=>({data:{id:x.label,label:x.one+x.number,side:l.includes(x)?'left':'right',rank:(l.includes(x)?l:r).indexOf(x),pep:(l.includes(x)?DATA.chains[p.left]:DATA.chains[p.right]).type==='peptide'}})),edges=es.map((e,i)=>({data:{id:'e'+i,source:e.a,target:e.b,type:cls(e.type),weight:e.fraction,label:(100*e.fraction).toFixed(1)+'%'}}));if(!cy){cy=cytoscape({container:$('#cy'),style:[{selector:'node',style:{label:'',color:'#07101e','font-size':11,'font-weight':800,'background-color':'#4bd7ff',width:38,height:38,'text-valign':'center','text-halign':'center'}},{selector:'node[pep]',style:{label:'data(label)','background-color':'#ffcf5b'}},{selector:'edge',style:{width:'mapData(weight, 0, 1, 2, 8)','line-color':'#94a3b8',opacity:.75,'curve-style':'bezier'}},{selector:'edge.bond',style:{'line-color':'#ff5578'}},{selector:'edge.ionic',style:{'line-color':'#c084fc'}},{selector:'edge.hydrophobic',style:{'line-color':'#7dff91'}},{selector:'.selected',style:{label:'data(label)','border-width':4,'border-color':'#fff',opacity:1}}]});cy.on('tap','node',e=>{let q=[...document.querySelectorAll('#sequence [data-r]')].find(x=>x.dataset.r===e.target.id());if(q)q.click()});cy.on('tap','edge',e=>{cy.elements().removeClass('selected');e.target.addClass('selected');e.target.connectedNodes().addClass('selected');$('#detail').innerHTML='<b>'+e.target.source().data('label')+' ↔ '+e.target.target().data('label')+'</b><br>'+e.target.data('label')+' geometric occupancy.'})}cy.elements().remove();cy.add(nodes);cy.add(edges);if(!nodes.length){cy.resize();return}let L=cy.nodes('[side = "left"]'),R=cy.nodes('[side = "right"]'),opts=(layout==='bipartite'||layout==='chain')?{name:'preset',positions:n=>({x:n.data('side')==='left'?180:1000,y:70+n.data('rank')*470/Math.max((n.data('side')==='left'?L:R).length-1,1)}),fit:true,padding:35}:{name:layout==='concentric'?'concentric':layout==='circle'?'circle':'cose',animate:false,fit:true,padding:35,concentric:n=>n.degree()};cy.layout(opts).run();requestAnimationFrame(()=>{cy.resize();cy.fit(cy.elements(),35)})}function refresh(){setTimeout(run,0)};['#pair','#occ','#type','#window','#start','#end'].forEach(s=>$(s).addEventListener('change',refresh));$('#occ').addEventListener('input',refresh);$('#layouts').querySelectorAll('button').forEach(b=>b.onclick=()=>{layout=b.dataset.layout;$('#layouts').querySelectorAll('button').forEach(x=>x.classList.toggle('on',x===b));run()});$('#tabs').addEventListener('click',e=>{if(e.target.dataset.v==='network')requestAnimationFrame(run)});if(window.ResizeObserver)new ResizeObserver(()=>{if(!$('#network').classList.contains('hidden'))refresh()}).observe($('#cy'));run()})();</script>'''


def build_interface_explorer(payload: str) -> str:
    """Add self-contained 3Dmol and synchronize the legacy view state.

    The older renderer is kept as a controlled fallback, but its selection
    function is patched before execution so all views share one residue key.
    """
    page = _build_interface_explorer_legacy(payload)
    threedmol_source = THREEDMOL_ASSET.read_text(encoding="utf-8")
    page = page.replace(
        "</script><script>const DATA=",
        "</script><script>" + threedmol_source + "</script><script>const DATA=",
        1,
    )
    # $3Dmol is already in the report; never request a portal-specific path.
    page = page.replace(
        "let s=document.createElement('script');s.src='/static/3Dmol-min.js';s.onload=load;s.onerror=load;document.head.appendChild(s)",
        "load()",
    )
    page = page.replace(
        "if(viewer)return;",
        "if(viewer){viewer.setStyle({},{cartoon:{color:'#3f5b76'}});let p=pair(),pep=[DATA.chains[p.left],DATA.chains[p.right]].find(c=>c.type==='peptide'),pepChain=p.left===pep?.label?p.left_chain:p.right_chain;viewer.setStyle({chain:pepChain},{cartoon:{color:'#ffcf5b'},stick:{radius:.18,color:'#ffdf85'}});let m=String(state.selected||'').split(':').at(-1).match(/([A-Z]{3})(-?\\d+)/);if(m)viewer.setStyle({resn:m[1],resi:m[2]},{stick:{radius:.24,color:'#ffffff'}});viewer.render();return;}",
        1,
    )
    page = page.replace(
        "viewer.addModel(atob(DATA.structure_pdb_b64),'pdb');viewer.setStyle",
        "viewer.addModel(atob(DATA.structure_pdb_b64),'pdb');viewer.setClickable({},{},atom=>{let p=pair(),candidate=[...DATA.chains[p.left].residues,...DATA.chains[p.right].residues].find(r=>r.number==atom.resi&&String(r.name||'').toUpperCase()===String(atom.resn||'').toUpperCase());if(candidate)window.dispatchEvent(new CustomEvent('pymacs-select',{detail:{residue:candidate.label}}));});viewer.setStyle",
        1,
    )
    bridge = """window.PyMACSSelection={residue:null,edge:null};window.addEventListener('pymacs-select',ev=>{const d=ev.detail||{};if(d.residue)select(d.residue,d.edge);});window.addEventListener('message',ev=>{const d=ev.data||{};if(d.type==='pymacs-rin-selection'&&d.residue)select(d.residue);});"""
    page = page.replace("function select(id,e){", bridge + "function select(id,e){", 1)
    page = page.replace(
        "state.selected=id;let ps=edges()",
        "state.selected=id;window.PyMACSSelection={residue:id,edge:e||null};document.dispatchEvent(new CustomEvent('pymacs-selection',{detail:window.PyMACSSelection}));let ps=edges()",
        1,
    )
    page = page.replace("window.parent?.postMessage({type:'pymacs-rin-selection',residue:id,edge:e&&[e.a,e.b]},'*');draw()}",
                        "window.parent?.postMessage({type:'pymacs-rin-selection',residue:id,edge:e&&[e.a,e.b]},'*');render()}", 1)
    # Rows are always peptide N→C irrespective of pair ordering.
    page = page.replace(
        "function contact(){let p=pair(),L=DATA.chains[p.left].residues,R=DATA.chains[p.right].residues,",
        "function contact(){let p=pair(),left=DATA.chains[p.left],right=DATA.chains[p.right],L=(left.type==='peptide'?left:right).residues,R=(left.type==='peptide'?right:left).residues,",
        1,
    )
    # Keep the wide receptor axis contained in its own scroller instead of
    # letting a large contact map crop the embedded dossier.
    page = page.replace(
        "$('#map').innerHTML='<p>Rows preserve peptide N→C order when a peptide is the left selected chain. Cell colour is geometric-contact occupancy.</p><table>",
        "$('#map').innerHTML='<p>Rows preserve peptide N→C order. Cell colour is geometric-contact occupancy.</p><div class=\"map-scroll\"><table>",
        1,
    )
    page = page.replace("+'</tr>').join('')+'</table>';$('#map')",
                        "+'</tr>').join('')+'</table></div>';$('#map')", 1)
    page = page.replace("#rin{display:none}",
                        "#map{overflow:hidden;max-height:610px}.map-scroll{max-width:100%;max-height:535px;overflow:auto;padding-bottom:6px;border:1px solid #29415b;border-radius:8px}.map-scroll table{min-width:max-content}.map-scroll th{position:sticky;top:0;background:#101d30;z-index:1}.map-scroll tr th:first-child{left:0;z-index:2}.map-scroll tr:first-child th:first-child{z-index:3}#rin{display:none}", 1)
    # SVG is the dependable primary network renderer in sandboxed dashboard
    # iframes.  Cytoscape remains vendored for standalone reports, but a
    # zero-sized GPU canvas must never hide the contact graph.
    page = page.replace(
        "#rin{display:none}",
        "#cy{display:none!important}#rin{display:block;background:#07111e;border-radius:10px}",
        1,
    )
    page = page.replace(
        ".res.sel{outline:3px solid var(--c)}",
        ".res.sel.pep{outline:3px solid #24e6ff;background:#083b4a;color:#efffff}.res.sel:not(.pep){outline:3px solid #ff5a78;background:#4b1726;color:#fff4f6}",
        1,
    )
    page = page.replace(
        "function select(id,e){",
        "function select(id,e){if(state.selected===id){state.selected=null;$('#detail').textContent='Selection cleared.';render();return}",
        1,
    )
    page = page.replace(
        "function draw(){let p=pair(),es=edges(),svg=$('#rin');",
        "function isPeptide(id){let p=pair();return[p.left,p.right].some(c=>DATA.chains[c].type==='peptide'&&DATA.chains[c].residues.some(r=>r.label===id))}function residue(id){let p=pair();return[...DATA.chains[p.left].residues,...DATA.chains[p.right].residues].find(r=>r.label===id)}function draw(){let p=pair(),es=edges(),svg=$('#rin');",
        1,
    )
    page = page.replace(
        "let A=[...new Set(es.map(e=>e.a))],B=[...new Set(es.map(e=>e.b))],pos={},at=(x,i,n)=>[x,80+i*470/Math.max(n-1,1)];A.forEach((x,i)=>pos[x]=at(180,i,A.length));B.forEach((x,i)=>pos[x]=at(1000,i,B.length));",
        "let A=[...new Set(es.map(e=>e.a))],B=[...new Set(es.map(e=>e.b))],all=[...A,...B],pos={},at=(x,i,n)=>[x,80+i*470/Math.max(n-1,1)];if(state.layout==='circle'||state.layout==='cose'){all.forEach((id,i)=>{let q=2*Math.PI*i/Math.max(all.length,1)-Math.PI/2;pos[id]=[600+330*Math.cos(q),305+230*Math.sin(q)]})}else if(state.layout==='concentric'){A.forEach((id,i)=>{let q=2*Math.PI*i/Math.max(A.length,1)-Math.PI/2;pos[id]=[600+360*Math.cos(q),305+240*Math.sin(q)]});B.forEach((id,i)=>{let q=2*Math.PI*i/Math.max(B.length,1)-Math.PI/2;pos[id]=[600+165*Math.cos(q),305+110*Math.sin(q)]})}else{A.forEach((x,i)=>pos[x]=at(180,i,A.length));B.forEach((x,i)=>pos[x]=at(1000,i,B.length))};let radius=Math.max(21,Math.min(38,220/Math.max(A.length,B.length,1)));",
        1,
    )
    page = page.replace(
        "r=\"${sel?32:25}\" fill=\"${sel?(isPeptide(id)?'#24e6ff':'#ff5a78'):(A.includes(id)?'#ffcf5b':'#4bd7ff')}\"/><text x=\"${q[0]}\" y=\"${q[1]+4}\" text-anchor=\"middle\" fill=\"#07101e\" font-weight=\"900\">${esc(id.split(':').at(-1))}</text>",
        "r=\"${sel?radius+5:radius}\" fill=\"${sel?(isPeptide(id)?'#24e6ff':'#ff5a78'):(A.includes(id)?'#ffcf5b':'#4bd7ff')}\"/><text x=\"${q[0]}\" y=\"${q[1]-12}\" text-anchor=\"middle\" fill=\"#07101e\" font-size=\"10\" font-weight=\"900\"><tspan x=\"${q[0]}\" dy=\"0\">${esc(String(residue(id)?.stable_id||'').split('|')[0]||'?')}</tspan><tspan x=\"${q[0]}\" dy=\"11\">${esc(residue(id)?.name||id.split(':').at(-1).replace(/-?\\d+$/,''))}</tspan><tspan x=\"${q[0]}\" dy=\"11\">${esc(residue(id)?.number||'')}</tspan></text>",
        1,
    )
    page = page.replace(
        "fill=\"${A.includes(id)?'#ffcf5b':'#4bd7ff'}\"",
        "fill=\"${sel?(isPeptide(id)?'#24e6ff':'#ff5a78'):(A.includes(id)?'#ffcf5b':'#4bd7ff')}\"",
        1,
    )
    page = page.replace("r=\"${sel?32:25}\"", "r=\"${sel?radius+5:radius}\"", 1)
    page = page.replace(
        "y=\"${q[1]+4}\" text-anchor=\"middle\" fill=\"#07101e\" font-weight=\"900\">${esc(id.split(':').at(-1))}</text>",
        "y=\"${q[1]-12}\" text-anchor=\"middle\" fill=\"#07101e\" font-size=\"10\" font-weight=\"900\"><tspan x=\"${q[0]}\" dy=\"0\">${esc(String(residue(id)?.stable_id||'').split('|')[0]||'?')}</tspan><tspan x=\"${q[0]}\" dy=\"11\">${esc(residue(id)?.name||id.split(':').at(-1).replace(/-?\\d+$/,''))}</tspan><tspan x=\"${q[0]}\" dy=\"11\">${esc(residue(id)?.number||'')}</tspan></text>",
        1,
    )
    page = page.replace(
        "c.fillStyle=isPep?'#ffcf5b':'#4bd7ff';",
        "c.fillStyle=state.selected===r.label?(isPep?'#24e6ff':'#ff5a78'):(isPep?'#ffcf5b':'#4bd7ff');",
        1,
    )
    page = page.replace("view:'network',a:.35}", "view:'network',a:.35,z:1,layout:'bipartite'}", 1)
    page = page.replace(
        ";$('#export').onclick=",
        ";$('#space').onwheel=e=>{e.preventDefault();state.z=Math.max(.35,Math.min(5,state.z*(e.deltaY<0?1.15:.87)));spatial()};$('#zoom-in').onclick=()=>{state.z=Math.min(5,state.z*1.3);spatial()};$('#zoom-out').onclick=()=>{state.z=Math.max(.35,state.z/1.3);spatial()};$('#zoom-reset').onclick=()=>{state.z=1;state.a=.35;spatial()};$('#layouts').querySelectorAll('button').forEach(b=>b.onclick=()=>{state.layout=b.dataset.layout;$('#layouts').querySelectorAll('button').forEach(x=>x.classList.toggle('on',x===b));draw()});$('#export').onclick=",
        1,
    )
    page = page.replace(
        "if(m)viewer.setStyle({resn:m[1],resi:m[2]},{stick:{radius:.24,color:'#ffffff'}});",
        "let selected=[...DATA.chains[p.left].residues,...DATA.chains[p.right].residues].find(r=>r.label===state.selected);if(selected){let colour=isPeptide(selected.label)?'#24e6ff':'#ff5a78';viewer.setStyle({resn:selected.name,resi:selected.number},{cartoon:{color:colour},stick:{radius:.32,color:colour},sphere:{scale:.38,color:colour}})};",
        1,
    )
    # The spatial view is contact-focused: all peptide residues remain for
    # sequence context, while protein labels are limited to active partners.
    page = page.replace(
        "function spatial(){let p=pair(),all=[...DATA.chains[p.left].residues,...DATA.chains[p.right].residues],by=new Map(all.map(r=>[r.label,r])),",
        "function spatial(){let p=pair(),raw=[...DATA.chains[p.left].residues,...DATA.chains[p.right].residues],involved=new Set(edges().flatMap(e=>[e.a,e.b])),all=raw.filter(r=>DATA.chains[p.left].type==='peptide'&&DATA.chains[p.left].residues.includes(r)||DATA.chains[p.right].type==='peptide'&&DATA.chains[p.right].residues.includes(r)||involved.has(r.label)),by=new Map(raw.map(r=>[r.label,r])),",
        1,
    )
    page = page.replace(
        "co=r=>r.position||[0,0,0],pr=r=>{let[x,y,z]=co(r),q=Math.cos(state.a),s=Math.sin(state.a);return[600+(x*q-z*s)*7,305-y*7+(x*s+z*q)]};",
        "co=r=>r.position||[0,0,0],bounds=all.map(co),centre=bounds.reduce((sum,v)=>sum.map((x,i)=>x+v[i]),[0,0,0]).map(x=>x/Math.max(bounds.length,1)),span=Math.max(1,...bounds.flatMap(v=>v.map((x,i)=>Math.abs(x-centre[i])))),scale=Math.min(11,250/span)*state.z,pr=r=>{let[x,y,z]=co(r);x-=centre[0];y-=centre[1];z-=centre[2];let q=Math.cos(state.a),s=Math.sin(state.a);return[600+(x*q-z*s)*scale,305-y*scale+(x*s+z*q)*scale*.16]};",
        1,
    )
    # Visually lift the peptide out of the receptor cartoon in the molecular
    # view.  The stable ID carries the PDB chain identifier.
    page = page.replace(
        "viewer.setStyle({},{cartoon:{color:'#53718c'}});viewer.zoomTo();viewer.render()",
        "let p=pair(),pep=[DATA.chains[p.left],DATA.chains[p.right]].find(c=>c.type==='peptide'),pepChain=p.left===pep?.label?p.left_chain:p.right_chain;viewer.setStyle({},{cartoon:{color:'#3f5b76'}});viewer.setStyle({chain:pepChain},{cartoon:{color:'#ffcf5b'},stick:{radius:.18,color:'#ffdf85'}});viewer.zoomTo();viewer.render()",
        1,
    )
    # Keep the molecular canvas bounded and retain a readable residue picker
    # beside it.  This avoids a 3Dmol canvas swallowing the whole tab in an
    # iframe and makes selections possible without scrolling away from the
    # structure.
    page = page.replace(
        "function structure(){let host=$('#structure');",
        "function structureInspector(){let p=pair(),left=DATA.chains[p.left],right=DATA.chains[p.right],pep=left.type==='peptide'?left:right,protein=left.type==='peptide'?right:left,partners=new Set(edges().flatMap(e=>[e.a,e.b])),proteinResidues=protein.residues.filter(r=>partners.has(r.label)),chips=(rs,pepChip)=>rs.map(r=>`<button class=\"res ${pepChip?'pep':''} ${state.selected===r.label?'sel':''}\" data-r=\"${esc(r.label)}\">${esc(r.one+r.number)}</button>`).join('');$('#structure-detail').innerHTML='<h3>Residue selector</h3><p>Peptide residues remain in N→C order. Protein entries are active contact partners only.</p><h4>Peptide · N→C</h4><div class=\"seq\">'+chips(pep.residues,true)+'</div><h4>Protein contact partners</h4><div class=\"seq\">'+(proteinResidues.length?chips(proteinResidues,false):'<span>No partners under the active filters.</span>')+'</div>';$('#structure-detail').querySelectorAll('[data-r]').forEach(n=>n.onclick=()=>select(n.dataset.r))}function structure(){let host=$('#mol');",
        1,
    )
    page = page.replace(
        "if(viewer){viewer.setStyle({},{cartoon:{color:'#3f5b76'}});let p=pair(),pep=[DATA.chains[p.left],DATA.chains[p.right]].find(c=>c.type==='peptide'),pepChain=p.left===pep?.label?p.left_chain:p.right_chain;viewer.setStyle({chain:pepChain},{cartoon:{color:'#ffcf5b'},stick:{radius:.18,color:'#ffdf85'}});let m=String(state.selected||'').split(':').at(-1).match(/([A-Z]{3})(-?\\d+)/);if(m)viewer.setStyle({resn:m[1],resi:m[2]},{stick:{radius:.24,color:'#ffffff'}});viewer.render();return;}",
        "if(viewer){viewer.setStyle({},{cartoon:{color:'#3f5b76'}});let p=pair(),pep=[DATA.chains[p.left],DATA.chains[p.right]].find(c=>c.type==='peptide'),pepChain=p.left===pep?.label?p.left_chain:p.right_chain;viewer.setStyle({chain:pepChain},{cartoon:{color:'#ffcf5b'},stick:{radius:.18,color:'#ffdf85'}});let m=String(state.selected||'').split(':').at(-1).match(/([A-Z]{3})(-?\\d+)/);if(m)viewer.setStyle({resn:m[1],resi:m[2]},{stick:{radius:.24,color:'#ffffff'}});requestAnimationFrame(()=>{viewer.resize();viewer.render()});return;}",
        1,
    )
    page = page.replace(
        "function render(){draw();sequence();contact();spatial();if(state.view==='structure')structure()}",
        "function render(){draw();sequence();contact();spatial();if(state.view==='structure'){structureInspector();structure()}}",
        1,
    )
    page = page.replace(
        "c.fillStyle='#fff';c.fillText(r.one+r.number,x+7,y-7)",
        "if(isPep||state.selected===r.label){c.fillStyle='#fff';c.fillText(r.one+r.number,x+7,y-7)}",
        1,
    )
    # A historical report can have no edge above its configured threshold.
    # Choose a bounded, occupancy-ranked initial threshold in that case; once
    # the user moves the slider, the active threshold is always respected.
    page = page.replace(
        "(()=>{let state={pair:DATA.pairs[0]?.id,occ:DATA.default_threshold||.5,",
        "(()=>{const CONTACT_EDGES=Array.isArray(DATA.edges)?DATA.edges:Object.values(DATA.edges||{}),defaultOcc=()=>{let requested=Number(DATA.default_threshold??.5),f=CONTACT_EDGES.map(e=>Number(e.fraction)).filter(Number.isFinite).sort((a,b)=>b-a);return !f.length||f.some(x=>x>=requested)?requested:f[Math.min(23,f.length-1)]};let state={pair:DATA.pairs[0]?.id,occ:defaultOcc(),",
        1,
    )
    page = page.replace("DATA.edges.forEach(e=>e.aggregate=e.fraction)", "CONTACT_EDGES.forEach(e=>e.aggregate=e.fraction)", 1)
    page = page.replace("return DATA.edges.filter(e=>p&&", "return CONTACT_EDGES.filter(e=>p&&", 1)
    # Windowed values are explicitly geometric; chemistry filters cannot claim
    # chemically validated occupancy without matching framewise evidence.
    page = page.replace("if(!times.length||!Object.keys(rows).length)", "if(!(DATA.presence?.available===true||(times.length&&Object.keys(rows).length)))", 1)
    page = page.replace("e.fraction=n?c/n:0});$('#note').textContent=", "e.fraction=n?c/n:0;e.window_geometric=$('#window').value!=='all'});$('#note').textContent=", 1)
    page = page.replace(
        "$('#window').onchange=e=>{$('#custom').classList.toggle('hidden',e.target.value!=='custom');recalc();render()};",
        "$('#window').onchange=e=>{$('#custom').classList.toggle('hidden',e.target.value!=='custom');let windowed=e.target.value!=='all';$('#type').disabled=windowed;if(windowed){state.kind='all';$('#type').value='all'}recalc();render()};",
        1,
    )
    page = page.replace("(state.kind==='all'||e.type===state.kind)", "(state.kind==='all'||(!e.window_geometric&&e.type===state.kind))")
    page = page.replace(
        "e.a.startsWith(p.left+':')&&e.b.startsWith(p.right+':')",
        "((e.a.startsWith(p.left+':')&&e.b.startsWith(p.right+':'))||(e.a.startsWith(p.right+':')&&e.b.startsWith(p.left+':')))",
        1,
    )
    return page


def build_cytoscape_interface_explorer(payload: str) -> str:
    """Build the report, then append the Cytoscape primary-renderer bridge."""
    page = build_interface_explorer(payload)
    # The contact-focused bipartite view is readable on first load.  The
    # exploratory force layout remains available as an explicit choice.
    page = page.replace('data-layout="cose" class="on"', 'data-layout="cose"', 1)
    page = page.replace('data-layout="bipartite">', 'data-layout="bipartite" class="on">', 1)
    return page.replace("</html>", CYTOSCAPE_BOOTSTRAP + "</html>")
