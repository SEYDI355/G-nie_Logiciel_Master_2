/* Fond 3D : un tardigrade (« ours d'eau ») qui tourne et marche.
   Réglages simples dans REGLAGES ci-dessous. Dépend de Three.js (chargé dans base.html). */
(function () {
  var canvas = document.getElementById('fond-3d');
  if (!canvas || !window.THREE) { return; }

  var REGLAGES = {
    couleur: 0x7fae7a,    // couleur du corps (vert mousse clair)
    vitesseRotation: 0.35, // tours du tardigrade (radians par seconde)
    vitesseMarche: 3.2,    // rythme des pattes
    taille: 1              // 1 = taille normale
  };

  var renderer;
  try {
    renderer = new THREE.WebGLRenderer({ canvas: canvas, alpha: true, antialias: true });
  } catch (e) { return; }                 // pas de WebGL : on laisse simplement le fond uni
  renderer.setClearColor(0x000000, 0);     // fond transparent
  renderer.setPixelRatio(Math.min(window.devicePixelRatio || 1, 2));

  var scene = new THREE.Scene();
  var camera = new THREE.PerspectiveCamera(35, 1, 0.1, 100);
  scene.add(new THREE.HemisphereLight(0xffffff, 0x8fb08c, 0.95));
  var soleil = new THREE.DirectionalLight(0xffffff, 0.85);
  soleil.position.set(3, 5, 4);
  scene.add(soleil);

  var matCorps = new THREE.MeshStandardMaterial({ color: REGLAGES.couleur, roughness: 0.5, metalness: 0 });
  var matGriffe = new THREE.MeshStandardMaterial({ color: 0xe6eee3, roughness: 0.6 });
  var matSombre = new THREE.MeshStandardMaterial({ color: 0x15231f, roughness: 0.4 });

  var tardi = new THREE.Group();
  scene.add(tardi);

  /* Corps : une suite de segments ellipsoïdaux qui se chevauchent */
  var segments = [];
  var N = 8;
  for (var i = 0; i < N; i++) {
    var t = i / (N - 1);
    var r = 0.55 + 0.45 * Math.sin(Math.PI * (0.15 + 0.7 * t));
    var seg = new THREE.Mesh(new THREE.SphereGeometry(1, 32, 24), matCorps);
    seg.scale.set(r * 0.95, r * 0.85, r);
    seg.position.set(-2.1 + 4.0 * t, 0, 0);
    tardi.add(seg);
    segments.push(seg);
  }

  /* Tête, bouche et yeux (le tardigrade regarde vers +x) */
  var tete = new THREE.Mesh(new THREE.SphereGeometry(1, 32, 24), matCorps);
  tete.scale.set(0.62, 0.55, 0.55);
  tete.position.set(2.4, 0.05, 0);
  tardi.add(tete);

  var bouche = new THREE.Mesh(new THREE.CylinderGeometry(0.13, 0.17, 0.25, 16), matSombre);
  bouche.rotation.z = Math.PI / 2;
  bouche.position.set(2.95, 0.03, 0);
  tardi.add(bouche);

  [-1, 1].forEach(function (cote) {
    var oeil = new THREE.Mesh(new THREE.SphereGeometry(0.07, 12, 10), matSombre);
    oeil.position.set(2.55, 0.3, cote * 0.3);
    tardi.add(oeil);
  });

  /* Quatre paires de pattes avec griffes */
  var pattes = [];
  [-1.4, -0.45, 0.5, 1.4].forEach(function (x, paire) {
    [-1, 1].forEach(function (cote) {
      var pivot = new THREE.Group();
      pivot.position.set(x, -0.5, cote * 0.62);
      pivot.rotation.x = cote * -0.55;                 // pattes inclinées vers l'extérieur

      var cuisse = new THREE.CylinderGeometry(0.2, 0.15, 0.7, 14);
      cuisse.translate(0, -0.35, 0);
      pivot.add(new THREE.Mesh(cuisse, matCorps));

      var pied = new THREE.Mesh(new THREE.SphereGeometry(0.19, 14, 12), matCorps);
      pied.scale.set(1, 0.8, 1.2);
      pied.position.set(0, -0.72, 0);
      pivot.add(pied);

      [-1, 0, 1].forEach(function (k) {
        var griffe = new THREE.Mesh(new THREE.ConeGeometry(0.04, 0.22, 8), matGriffe);
        griffe.position.set(0.08, -0.86, k * 0.09);
        griffe.rotation.z = Math.PI + 0.35;            // pointe vers le bas et l'avant
        pivot.add(griffe);
      });

      tardi.add(pivot);
      // Marche : les pattes opposées sont en décalage de phase
      pattes.push({ pivot: pivot, phase: paire * 0.9 + (cote > 0 ? Math.PI : 0) });
    });
  });

  tardi.rotation.z = 0.12;                              // léger penché

  function ajuster() {
    var l = window.innerWidth, h = window.innerHeight, ratio = l / h;
    renderer.setSize(l, h, false);
    camera.aspect = ratio;
    camera.position.set(0, 1.4, ratio < 1 ? 9 / ratio : 9);
    camera.lookAt(0, 0, 0);
    camera.updateProjectionMatrix();
    tardi.scale.setScalar(REGLAGES.taille);
    tardi.position.x = ratio > 1.2 ? -0.8 : 0;
    renderer.render(scene, camera);
  }
  window.addEventListener('resize', ajuster);
  ajuster();

  /* Respect de « réduire les animations » : image fixe */
  if (window.matchMedia && window.matchMedia('(prefers-reduced-motion: reduce)').matches) { return; }

  var horloge = new THREE.Clock();
  function boucle() {
    var dt = horloge.getDelta();
    var t = horloge.elapsedTime;

    tardi.rotation.y += dt * REGLAGES.vitesseRotation;  // il tourne sur lui-même
    tardi.position.y = Math.sin(t * 1.3) * 0.12;        // petit balancement vertical

    segments.forEach(function (s, i) {                  // ondulation du corps
      s.position.y = Math.sin(t * 2.2 - i * 0.7) * 0.05;
    });
    pattes.forEach(function (p) {                       // il marche
      p.pivot.rotation.z = Math.sin(t * REGLAGES.vitesseMarche + p.phase) * 0.45;
    });

    renderer.render(scene, camera);
    requestAnimationFrame(boucle);
  }
  boucle();
})();
