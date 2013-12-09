Map =
  center: ->
    x = Math.floor(world.realSize.x / 2)
    y = Math.floor(world.realSize.y / 2)
    game.camera.unconstrain().focus(x, y)