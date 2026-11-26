/* Intégration Cytoscape.js pour l'explorateur métabolique KEGG (pmeo) */
(function () {
  var tabEl = document.getElementById("tab-network");
  if (!tabEl || typeof window.cytoscape !== "function") {
    return;
  }

  var keggUrl = tabEl.dataset.keggUrl;
  var measurementsUrl = tabEl.dataset.measurementsUrl;
  var cy = null;
  var rawKegg = null;
  var datasetPayload = null;
  var initialized = false;
  var loading = false;

  function el(id) {
    return document.getElementById(id);
  }

  function escapeHtml(text) {
    return String(text == null ? "" : text)
      .replace(/&/g, "&amp;")
      .replace(/</g, "&lt;")
      .replace(/>/g, "&gt;")
      .replace(/"/g, "&quot;")
      .replace(/'/g, "&#39;");
  }

  function fmtNum(val, digits) {
    if (val === null || val === undefined || isNaN(val)) {
      return "—";
    }
    return Number(val).toFixed(digits === undefined ? 2 : digits);
  }

  function fmtSci(val) {
    if (val === null || val === undefined || isNaN(val)) {
      return "—";
    }
    var n = Number(val);
    return n < 0.001 && n > 0 ? n.toExponential(2) : n.toFixed(4);
  }

  function lerpColor(hexA, hexB, t) {
    t = Math.max(0, Math.min(1, t));
    var a = parseInt(hexA.slice(1), 16);
    var b = parseInt(hexB.slice(1), 16);
    var r = Math.round(((a >> 16) & 255) + (((b >> 16) & 255) - ((a >> 16) & 255)) * t);
    var g = Math.round(((a >> 8) & 255) + (((b >> 8) & 255) - ((a >> 8) & 255)) * t);
    var bl = Math.round((a & 255) + ((b & 255) - (a & 255)) * t);
    return "rgb(" + r + "," + g + "," + bl + ")";
  }

  function seqColor(val, minV, maxV) {
    var t = maxV === minV ? 0.65 : (val - minV) / (maxV - minV);
    return t < 0.5
      ? lerpColor("#fef08a", "#f97316", t * 2)
      : lerpColor("#f97316", "#991b1b", (t - 0.5) * 2);
  }

  function divColor(val, maxAbs) {
    var lim = maxAbs > 0 ? maxAbs : 1;
    var clamped = Math.max(-lim, Math.min(lim, val));
    if (clamped >= 0) {
      return lerpColor("#f8fafc", "#dc2626", clamped / lim);
    }
    return lerpColor("#f8fafc", "#2563eb", Math.abs(clamped) / lim);
  }

  function createCytoscapeInstance() {
    cy = window.cytoscape({
      container: el("network-cy"),
      style: [
        {
          selector: 'node[node_type = "metabolite"]',
          style: {
            "background-color": "data(viz_color)",
            "width": "data(viz_size)",
            "height": "data(viz_size)",
            "border-width": "data(viz_border_w)",
            "border-color": "data(viz_border_c)",
            "label": "data(display_label)",
            "font-size": "10px",
            "color": "#15231f",
            "text-valign": "bottom",
            "text-margin-y": 4,
            "text-background-color": "#ffffff",
            "text-background-opacity": 0.8,
            "text-background-padding": "2px"
          }
        },
        {
          selector: 'node[node_type = "reaction"]',
          style: {
            "shape": "round-rectangle",
            "background-color": "#94a3b8",
            "label": "data(label)",
            "font-size": "9px",
            "color": "#ffffff",
            "text-valign": "center",
            "text-halign": "center",
            "width": 48,
            "height": 20
          }
        },
        {
          selector: 'node[evidence_source = "KEGG_pmeo"]',
          style: { "background-color": "#16a34a" }
        },
        {
          selector: "edge",
          style: {
            "width": 1.6,
            "line-color": "#94a3b8",
            "target-arrow-color": "#94a3b8",
            "target-arrow-shape": "triangle",
            "curve-style": "bezier",
            "arrow-scale": 0.85
          }
        },
        {
          selector: 'edge[evidence_source = "KEGG_pmeo"]',
          style: {
            "line-color": "#16a34a",
            "target-arrow-color": "#16a34a",
            "source-arrow-color": "#16a34a"
          }
        },
        {
          selector: "edge[?reversible]",
          style: { "source-arrow-shape": "triangle", "source-arrow-color": "#94a3b8" }
        },
        {
          selector: ".highlighted",
          style: {
            "border-width": 3.5,
            "border-color": "#2563eb",
            "line-color": "#2563eb",
            "target-arrow-color": "#2563eb",
            "source-arrow-color": "#2563eb",
            "z-index": 999
          }
        },
        {
          selector: ".dimmed",
          style: { "opacity": 0.16 }
        }
      ]
    });

    bindCyEvents();
  }

  function applyVisualMapping() {
    if (!cy) return;
    var colorMode = el("netColorMode").value;
    var sizeMode = el("netSizeMode").value;
    var measurements = (datasetPayload && datasetPayload.measurements) || {};
    var metNodes = cy.nodes('[node_type = "metabolite"]');

    var vals = [];
    var fcs = [];
    var degs = [];
    var betws = [];

    metNodes.forEach(function (n) {
      var m = measurements[n.id()];
      if (m) {
        if (m.value !== null && m.value !== undefined) vals.push(Number(m.value));
        if (m.log2fc !== null && m.log2fc !== undefined) fcs.push(Number(m.log2fc));
      }
      degs.push(Number(n.data("degree") || 1));
      betws.push(Number(n.data("betweenness") || 0));
    });

    var minVal = vals.length ? Math.min.apply(null, vals) : 0;
    var maxVal = vals.length ? Math.max.apply(null, vals) : 1;
    var maxAbsFc = fcs.length ? Math.max.apply(null, fcs.map(Math.abs)) : 2;
    var maxDeg = degs.length ? Math.max.apply(null, degs) : 1;
    var maxBetw = betws.length ? Math.max.apply(null, betws) : 0.01;

    cy.batch(function () {
      metNodes.forEach(function (n) {
        var m = measurements[n.id()] || null;
        var v = m ? m.value : null;
        var fc = m ? m.log2fc : null;
        var isMeasured = m !== null && ((v !== null && v !== undefined) || (fc !== null && fc !== undefined));

        n.data("measurement", m);
        n.data("quant_value", v);
        n.data("log2fc", fc);
        n.data("evidence_status", isMeasured ? "measured" : "inferred_kegg");

        var baseLabel = n.data("label");
        if (m && m.metabolite_name) {
          n.data("display_label", m.metabolite_name);
        } else if (baseLabel && baseLabel !== "nan") {
          n.data("display_label", baseLabel);
        } else {
          n.data("display_label", n.id());
        }

        var color = isMeasured ? "#f97316" : "#cbd5e1";
        if (colorMode === "value") {
          color = v !== null && v !== undefined ? seqColor(Number(v), minVal, maxVal) : "#e2e8f0";
        } else if (colorMode === "log2fc") {
          color = fc !== null && fc !== undefined ? divColor(Number(fc), maxAbsFc) : "#e2e8f0";
        } else if (colorMode === "betweenness") {
          color = seqColor(Number(n.data("betweenness") || 0), 0, maxBetw || 1);
        }

        var size = 22;
        if (sizeMode === "value" && v !== null && v !== undefined) {
          var tv = maxVal === minVal ? 0.6 : (Number(v) - minVal) / (maxVal - minVal);
          size = 18 + tv * 30;
        } else if (sizeMode === "log2fc" && fc !== null && fc !== undefined) {
          var tf = maxAbsFc > 0 ? Math.abs(Number(fc)) / maxAbsFc : 0.5;
          size = 18 + tf * 30;
        } else if (sizeMode === "degree") {
          var td = Number(n.data("degree") || 1) / (maxDeg || 1);
          size = 18 + td * 28;
        } else if (sizeMode !== "fixed") {
          size = 15;
        }

        n.data("viz_color", color);
        n.data("viz_size", Math.round(size));
        n.data("viz_border_w", isMeasured ? 2.5 : 1);
        n.data("viz_border_c", isMeasured ? "#15231f" : "#64748b");
      });
    });

    updateLegend(colorMode, minVal, maxVal, maxAbsFc, maxBetw, vals.length, fcs.length);
  }

  function updateLegend(colorMode, minVal, maxVal, maxAbsFc, maxBetw, nVal, nFc) {
    var box = el("net-legend");
    if (colorMode === "value") {
      box.innerHTML =
        "<span><strong>Abondance (" + nVal + " mesurés)</strong>" +
        '<div class="net-grad-bar" style="background:linear-gradient(to right,#fef08a,#f97316,#991b1b);"></div>' +
        '<div class="net-grad-labels"><span>' + minVal.toFixed(0) + "</span><span>" + maxVal.toFixed(0) + "</span></div></span>" +
        '<span><i class="swatch" style="background:#e2e8f0;border:1px solid #64748b;"></i>Inféré KEGG (non mesuré)</span>';
    } else if (colorMode === "log2fc") {
      box.innerHTML =
        "<span><strong>Log2 Fold Change (" + nFc + " mesurés)</strong>" +
        '<div class="net-grad-bar" style="background:linear-gradient(to right,#2563eb,#f8fafc,#dc2626);"></div>' +
        '<div class="net-grad-labels"><span>-' + maxAbsFc.toFixed(1) + "</span><span>0</span><span>+" + maxAbsFc.toFixed(1) + "</span></div></span>" +
        '<span><i class="swatch" style="background:#e2e8f0;border:1px solid #64748b;"></i>Non mesuré</span>' +
        '<span><i class="swatch" style="background:#16a34a;border-radius:3px;"></i>Réaction pmeo</span>';
    } else if (colorMode === "betweenness") {
      box.innerHTML =
        "<span><strong>Centralité NetworkX</strong>" +
        '<div class="net-grad-bar" style="background:linear-gradient(to right,#fef08a,#f97316,#991b1b);"></div>' +
        '<div class="net-grad-labels"><span>0</span><span>' + maxBetw.toFixed(4) + "</span></div></span>";
    } else {
      box.innerHTML =
        '<span><i class="swatch" style="background:#f97316;border:2px solid #15231f;"></i>Métabolite mesuré</span>' +
        '<span><i class="swatch" style="background:#cbd5e1;border:1px solid #64748b;"></i>Métabolite inféré (KEGG)</span>' +
        '<span><i class="swatch" style="background:#16a34a;border-radius:3px;"></i>Réaction attestée (KEGG pmeo)</span>' +
        '<span><i class="swatch" style="background:#94a3b8;border-radius:3px;"></i>Réaction de référence</span>';
    }
  }

  function updateSummaryText(nMet, nRxn, nEdges) {
    var summary = (datasetPayload && datasetPayload.summary) || {};
    var parts = [
      "Sous-graphe affiché : " + nMet + " métabolites, " + nRxn + " réactions, " + nEdges + " arêtes"
    ];
    if (summary.total_rows !== undefined) {
      parts.push(
        "Jeu de données (" + (summary.selected_condition || "—") +
        " vs " + (summary.control_condition || "—") + ") : " +
        summary.mapped_count + "/" + summary.total_rows + " mappés KEGG, " +
        summary.in_graph_count + " présents dans le réseau pmeo, " +
        summary.significant_count + " significatifs"
      );
    }
    el("net-summary").textContent = parts.join(" · ");
  }

  function renderGraph() {
    if (!rawKegg || !cy) return;

    var mode = el("netGraphMode").value;
    var selectedPw = el("netPathwaySelect").value;
    var hideCurrency = el("netHideCurrency").checked;
    var onlyPmeo = el("netOnlyPmeo").checked;
    var layoutName = el("netLayoutSelect").value;

    var sourceSet = rawKegg[mode];
    var currencyIds = new Set(
      sourceSet.nodes.filter(function (n) { return n.data.is_currency; }).map(function (n) { return n.data.id; })
    );

    var edges = sourceSet.edges.filter(function (e) {
      var inPw = selectedPw === "ALL" || (e.data.pathways || []).indexOf(selectedPw) !== -1;
      if (!inPw) return false;
      if (hideCurrency && (currencyIds.has(e.data.source) || currencyIds.has(e.data.target))) return false;
      if (onlyPmeo && mode === "monopartite" && e.data.evidence_source !== "KEGG_pmeo") return false;
      return true;
    });

    if (onlyPmeo && mode === "bipartite") {
      var pmeoRxnIds = new Set(
        sourceSet.nodes
          .filter(function (n) { return n.data.node_type === "reaction" && n.data.evidence_source === "KEGG_pmeo"; })
          .map(function (n) { return n.data.id; })
      );
      edges = edges.filter(function (e) {
        return pmeoRxnIds.has(e.data.source) || pmeoRxnIds.has(e.data.target);
      });
    }

    var activeIds = new Set();
    edges.forEach(function (e) {
      activeIds.add(e.data.source);
      activeIds.add(e.data.target);
    });

    var hasPresetCoords = false;
    var nodes = sourceSet.nodes
      .filter(function (n) { return activeIds.has(n.data.id); })
      .map(function (n) {
        var clone = JSON.parse(JSON.stringify(n));
        var posMap = clone.data.pos_by_pathway || {};
        if (selectedPw !== "ALL" && posMap[selectedPw]) {
          clone.position = { x: posMap[selectedPw].x, y: posMap[selectedPw].y };
          hasPresetCoords = true;
        }
        clone.data.viz_color = "#cbd5e1";
        clone.data.viz_size = 22;
        clone.data.viz_border_w = 1;
        clone.data.viz_border_c = "#64748b";
        clone.data.display_label = clone.data.label && clone.data.label !== "nan" ? clone.data.label : clone.data.id;
        return clone;
      });

    var clonedEdges = edges.map(function (e) { return JSON.parse(JSON.stringify(e)); });

    cy.elements().remove();
    cy.add(nodes);
    cy.add(clonedEdges);

    applyVisualMapping();

    if (layoutName === "preset" && !hasPresetCoords) {
      layoutName = "cose";
    }

    var layoutOpts = { name: layoutName, animate: false, padding: 30 };
    if (layoutName === "cose") {
      layoutOpts = { name: "cose", animate: false, nodeRepulsion: 9000, idealEdgeLength: 75 };
    } else if (layoutName === "concentric") {
      layoutOpts = {
        name: "concentric",
        animate: false,
        concentric: function (n) { return (n.data("betweenness") || 0) * 1000 + (n.data("degree") || 0); },
        levelWidth: function () { return 2; }
      };
    }

    cy.resize();
    cy.layout(layoutOpts).run();
    cy.fit(undefined, 35);

    var nMet = cy.nodes('[node_type = "metabolite"]').length;
    var nRxn = cy.nodes('[node_type = "reaction"]').length;
    updateSummaryText(nMet, nRxn, clonedEdges.length);
  }

  function bindCyEvents() {
    cy.on("tap", "node", function (evt) {
      var node = evt.target;
      var d = node.data();

      cy.elements().removeClass("highlighted dimmed");
      var hood = node.closedNeighborhood();
      cy.elements().difference(hood).addClass("dimmed");
      hood.addClass("highlighted");

      if (d.node_type === "metabolite") {
        var m = d.measurement;
        var badgeColor = m ? "#3e7b45" : "#64748b";
        var html =
          '<span class="net-badge" style="background:' + badgeColor + ';">Métabolite (' +
          escapeHtml(m ? "mesuré : " + m.significance : "inféré KEGG") + ")</span>" +
          '<div class="net-kv"><strong>ID KEGG</strong><a href="https://www.kegg.jp/entry/' +
          encodeURIComponent(d.id) + '" target="_blank" rel="noopener">' + escapeHtml(d.id) + " ↗</a></div>" +
          '<div class="net-kv"><strong>Nom</strong>' + escapeHtml((m && m.metabolite_name) || d.full_name || d.id) + "</div>";

        if (m) {
          html +=
            '<div class="net-kv"><strong>Mesures (' + escapeHtml(m.condition) + " vs " + escapeHtml(m.control_condition) + ")</strong>" +
            "Abondance : <b>" + fmtNum(m.value, 1) + "</b> (ctrl : " + fmtNum(m.control_value, 1) + ")<br>" +
            "Log2FC : <b>" + fmtNum(m.log2fc, 3) + "</b><br>" +
            "p-value : " + fmtSci(m.pvalue) + " · BH : " + fmtSci(m.padj_bh) + "<br>" +
            "Mapping : " + escapeHtml(m.mapping_status) + (m.pathway_label ? " (" + escapeHtml(m.pathway_label) + ")" : "") +
            "</div>";
        }
        html +=
          '<div class="net-kv"><strong>Topologie NetworkX</strong>Degré : ' + escapeHtml(d.degree) +
          " (in: " + escapeHtml(d.in_degree) + ", out: " + escapeHtml(d.out_degree) + ") · Betweenness : " + escapeHtml(d.betweenness) + "</div>" +
          '<div class="net-kv"><strong>Voies KEGG (' + (d.pathways || []).length + ")</strong>" + escapeHtml((d.pathways || []).join(", ")) + "</div>";

        el("net-details").innerHTML = html;
      } else {
        var rxnColor = d.evidence_source === "KEGG_pmeo" ? "#16a34a" : "#64748b";
        el("net-details").innerHTML =
          '<span class="net-badge" style="background:' + rxnColor + ';">Réaction (' + escapeHtml(d.evidence_source) + ")</span>" +
          '<div class="net-kv"><strong>ID KEGG</strong><a href="https://www.kegg.jp/entry/' +
          encodeURIComponent(d.id) + '" target="_blank" rel="noopener">' + escapeHtml(d.id) + " ↗</a></div>" +
          '<div class="net-kv"><strong>Gènes pmeo (' + escapeHtml(d.gene_count || 0) + ")</strong>" +
          escapeHtml((d.genes_pmeo || []).join(", ") || "Aucun détecté") + "</div>" +
          '<div class="net-kv"><strong>KO / EC</strong>KO : ' + escapeHtml((d.ko_ids || []).join(", ") || "—") +
          "<br>EC : " + escapeHtml((d.ec_ids || []).join(", ") || "—") + "</div>" +
          '<div class="net-kv"><strong>Réversible</strong>' + (d.reversible ? "Oui" : "Non") + "</div>" +
          '<div class="net-kv"><strong>Voies KEGG</strong>' + escapeHtml((d.pathways || []).join(", ")) + "</div>";
      }
    });

    cy.on("tap", "edge", function (evt) {
      var d = evt.target.data();
      el("net-details").innerHTML =
        '<span class="net-badge" style="background:#334155;">Arête (' + escapeHtml(d.role || "projection monopartite") + ")</span>" +
        '<div class="net-kv"><strong>Transition</strong>' + escapeHtml(d.source) + " → " + escapeHtml(d.target) + "</div>" +
        '<div class="net-kv"><strong>Réactions associées</strong>' + escapeHtml((d.reactions || []).join(", ") || "Lien direct biparti") + "</div>" +
        '<div class="net-kv"><strong>Gènes pmeo</strong>' + escapeHtml((d.genes_pmeo || []).join(", ") || "—") + "</div>" +
        '<div class="net-kv"><strong>Réversible</strong>' + (d.reversible ? "Oui" : "Non") + "</div>";
    });

    cy.on("tap", function (evt) {
      if (evt.target === cy) {
        cy.elements().removeClass("highlighted dimmed");
      }
    });
  }

  function populateControls() {
    var pwSel = el("netPathwaySelect");
    pwSel.innerHTML = "";
    (rawKegg.pathways || []).forEach(function (pw) {
      var opt = document.createElement("option");
      opt.value = pw.id;
      opt.textContent = pw.id + " — " + pw.title;
      pwSel.appendChild(opt);
    });
    var allOpt = document.createElement("option");
    allOpt.value = "ALL";
    allOpt.textContent = "Tout le réseau (lourd)";
    pwSel.appendChild(allOpt);
    if ((rawKegg.pathways || []).some(function (p) { return p.id === "pmeo00500"; })) {
      pwSel.value = "pmeo00500";
    }

    var condSel = el("netConditionSelect");
    condSel.innerHTML = "";
    var conditions = (datasetPayload && datasetPayload.summary && datasetPayload.summary.conditions) || [];
    if (!conditions.length) {
      var noneOpt = document.createElement("option");
      noneOpt.value = "";
      noneOpt.textContent = "Aucune comparaison";
      condSel.appendChild(noneOpt);
      el("netColorMode").value = "default";
      el("netSizeMode").value = "degree";
    } else {
      conditions.forEach(function (c) {
        var opt = document.createElement("option");
        opt.value = c;
        opt.textContent = c;
        condSel.appendChild(opt);
      });
      if (datasetPayload.summary.selected_condition) {
        condSel.value = datasetPayload.summary.selected_condition;
      }
    }
  }

  function bindUiControls() {
    ["netGraphMode", "netPathwaySelect", "netLayoutSelect", "netHideCurrency", "netOnlyPmeo"].forEach(function (id) {
      el(id).addEventListener("change", renderGraph);
    });
    ["netColorMode", "netSizeMode"].forEach(function (id) {
      el(id).addEventListener("change", applyVisualMapping);
    });
    el("netConditionSelect").addEventListener("change", function (e) {
      var cond = e.target.value;
      if (!cond) return;
      fetch(measurementsUrl + "?condition=" + encodeURIComponent(cond))
        .then(function (r) { return r.json(); })
        .then(function (data) {
          datasetPayload = data;
          applyVisualMapping();
          var nMet = cy.nodes('[node_type = "metabolite"]').length;
          var nRxn = cy.nodes('[node_type = "reaction"]').length;
          var nEdges = cy.edges().length;
          updateSummaryText(nMet, nRxn, nEdges);
        });
    });
    el("netFitBtn").addEventListener("click", function () {
      if (cy) cy.fit(undefined, 35);
    });
    el("netExportBtn").addEventListener("click", function () {
      if (!cy) return;
      var a = document.createElement("a");
      a.href = cy.png({ full: true, scale: 2, bg: "#ffffff" });
      a.download = "pmeo_" + el("netPathwaySelect").value + ".png";
      a.click();
    });
    el("netSearchInput").addEventListener("input", function (e) {
      if (!cy) return;
      var q = e.target.value.trim().toLowerCase();
      cy.elements().removeClass("highlighted dimmed");
      if (!q) return;
      var matches = cy.nodes().filter(function (n) {
        var m = n.data("measurement");
        var blob = [
          n.id(),
          n.data("full_name") || "",
          n.data("display_label") || "",
          m ? m.metabolite_name : "",
          (n.data("genes_pmeo") || []).join(" "),
          (n.data("ko_ids") || []).join(" ")
        ].join(" ").toLowerCase();
        return blob.indexOf(q) !== -1;
      });
      if (matches.length) {
        cy.elements().difference(matches.closedNeighborhood()).addClass("dimmed");
        matches.addClass("highlighted");
        cy.animate({ fit: { eles: matches, padding: 80 } }, { duration: 250 });
      }
    });
  }

  function initNetwork() {
    if (initialized || loading) return;
    loading = true;
    el("net-summary").textContent = "Chargement du réseau métabolique KEGG et des mesures...";

    Promise.all([
      fetch(keggUrl).then(function (r) {
        if (!r.ok) throw new Error("Erreur HTTP KEGG (" + r.status + ")");
        return r.json();
      }),
      fetch(measurementsUrl).then(function (r) {
        if (!r.ok) throw new Error("Erreur HTTP mesures (" + r.status + ")");
        return r.json();
      })
    ])
      .then(function (results) {
        rawKegg = results[0];
        datasetPayload = results[1];
        createCytoscapeInstance();
        populateControls();
        bindUiControls();
        initialized = true;
        loading = false;
        renderGraph();
      })
      .catch(function (err) {
        loading = false;
        el("net-summary").textContent = "Impossible de charger le réseau : " + err.message;
      });
  }

  window.TardigradeNetwork = {
    onTabShown: function () {
      if (!initialized) {
        initNetwork();
      } else if (cy) {
        cy.resize();
        cy.fit(undefined, 35);
      }
    }
  };
})();
