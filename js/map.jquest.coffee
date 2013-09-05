Map =
  center: ->
    x = game.renderer.width / 2
    y = game.renderer.height / 2
    width = (world.size.x * 32) / 2
    height = (world.size.y * 32) / 2
    newx = (x - width) - world.position.x
    newy = (y - height) - world.position.y
    world.pan(newx, newy)