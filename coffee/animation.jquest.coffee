Animation =

  loadChar: (name, eight = true) ->

    textures = []
    directions = ['s','w','e','n','sw','se','nw','ne']
    baseTexture = game.cache.getTexture(name).baseTexture
    sprite = new gf.Sprite()
    width = 32
    height = 32

    for y in [0...8]
      textures[directions[y]] = []
      for x in [0...6]
        pointX = x * 32
        pointY = y * 32
        textures[directions[y]].push new gf.Texture(baseTexture, new gf.Rectangle(pointX, pointY, width, height))

    (sprite.addAnimation(k, v, 0.08, true) for k, v of textures)
    sprite.direction = false
    sprite


