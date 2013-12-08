$(window).ready(function() {
  var engine, tileset;
  engine = new Joy.Engine(document.getElementById('game'));
  tileset = new Joy.Tileset({
    src: "resources/img/free_tileset_version_10.png",
    width: 32,
    height: 32
  });
  return $.getJSON('resources/map/island2.json', function(json) {
    return engine.createScene(function(scene) {
      var layer, _i, _len, _ref;
      _ref = json.layers;
      for (_i = 0, _len = _ref.length; _i < _len; _i++) {
        layer = _ref[_i];
        layer['lines'] = layer.width;
        layer['columns'] = layer.height;
        layer['tileset'] = tileset;
        scene.addChild(new Joy.Tilemap(layer));
      }
      return false;
    });
  });
});
