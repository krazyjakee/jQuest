Animation =

  loadChar: (name, eight = true) ->

    textures = []
    directions = ['s','w','e','n','sw','se','nw','ne']
    bt = gf.BaseTexture 'resources/img/'+name+'.png'
    width = 32
    height = 32

    for y in [0...8]
      for x in [0...6]
        pointX = x * 32
        pointY = y * 32
        textures[directions[y]].push new gf.Texture(bt, gf.Rectangle(pointX, pointY, width, height))

    return gf.Sprite
      n:
        frames: textures['n']
        rate: 2
      s:
        fames: textures['s']
        rate: 2
      w:
        fames: textures['w']
        rate: 2
      e:
        fames: textures['e']
        rate: 2
      nw:
        fames: textures['nw']
        rate: 2
      ne:
        fames: textures['ne']
        rate: 2
      sw:
        fames: textures['sw']
        rate: 2
      se:
        fames: textures['se']
        rate: 2
    , 1000


