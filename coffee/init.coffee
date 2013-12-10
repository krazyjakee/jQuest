world = game = false

$(window).ready ->
    
  gf.Texture.SCALE_MODE.DEFAULT = gf.Texture.SCALE_MODE.NEAREST

  game = new gf.Game 'game',
    width: window.innerWidth
    height: window.innerHeight
    background: 0x000000
    renderer: gf.RENDERER.AUTO

  game.load.tilemap 'island2', 'resources/map/island2.json', null, gf.FILE_FORMAT.JSON

  game.load.once 'complete', ->
    game.render()

  game.load.start()

window.onresize = ->
  #game.resize(window.innerWidth, window.innerHeight)
  #Map.center()