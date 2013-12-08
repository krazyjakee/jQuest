$(window).ready ->
  engine = new Joy.Engine document.getElementById('game')
  
  tileset = new Joy.Tileset
    src: "resources/img/free_tileset_version_10.png"
    width: 32
    height: 32

  $.getJSON 'resources/map/island2.json', (json) ->
    engine.createScene (scene) -> 
      for layer in json.layers
        layer['lines'] = layer.width
        layer['columns'] = layer.height
        layer['tileset'] = tileset
        scene.addChild(new Joy.Tilemap(layer))
      false