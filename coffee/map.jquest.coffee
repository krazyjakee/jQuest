Map = 

  tileClick: (e) ->
    tile = e.target
    tile.hitArea = new gf.Rectangle(-32, -32, 64, 64)
    tile.mass = Infinity
    tile.enablePhysics game.physics
    tile.alpha = 0.5
    console.log tile