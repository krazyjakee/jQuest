game = false

$(window).ready ->
    
  gf.Texture.SCALE_MODE.DEFAULT = gf.Texture.SCALE_MODE.NEAREST

  game = new gf.Game 'game',
    width: 1024
    height: 640
    background: 0x000000
    renderer: gf.RENDERER.AUTO

  game.load.tilemap 'island2', 'resources/map/island2.json', null, gf.FILE_FORMAT.JSON
  game.load.image '6Actor_5', 'resources/img/6Actor_5.png', null, gf.ATLAS_FORMAT.JSON_HASH

  game.load.once 'complete', ->
    playerChar = Animation.loadChar('6Actor_5')
    game.world.add.tilemap('island2', true)
    game.world.addChild(playerChar)
    playerChar.goto(0, 's').play()
    game.render()

  game.load.start()