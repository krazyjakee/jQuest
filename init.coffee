world = game = false

$(window).ready ->
    
  game = new gf.Game 'game',
    width: window.innerWidth
    height: window.innerHeight
    background: 0x000000
    interactive: true

  game.loader.on 'complete', ->
    game.loadWorld 'island2'
    world = game.world
    Map.center()
    game.render()

  game.loader.on 'progress', (e) ->
    false

  game.loader.load [
    name: 'island2'
    src: 'resources/map/island2.json'
  ]

window.onresize = ->
  game.renderer.resize(window.innerWidth, window.innerHeight)
  Map.center()