Animation =

  loadChar: (name, eight = true) ->

    textures = []
    directions = ['s','w','e','n','sw','se','nw','ne']
    baseTexture = game.cache.getTexture(name).baseTexture
    width = 32
    height = 32

    for y in [0...8]
      textures[directions[y]] = []
      for x in [0...6]
        pointX = x * 32
        pointY = y * 32
        textures[directions[y]].push new gf.Texture(baseTexture, gf.Rectangle(pointX, pointY, width, height))

    return gf.Sprite
      n:
        frames: textures['n']
      s:
        frames: textures['s']
      w:
        frames: textures['w']
      e:
        frames: textures['e']
      nw:
        frames: textures['nw']
      ne:
        frames: textures['ne']
      sw:
        frames: textures['sw']
      se:
        frames: textures['se']
    , 1000, 's'


