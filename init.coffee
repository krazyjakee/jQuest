world = game = false

$(window).ready ->
    
  game = new gf.Game 'game',
    width: window.innerWidth
    height: window.innerHeight
    background: 0x000000
    interactive: true

  game.on 'complete', ->
    game.loadWorld 'island2'
    world = game.world
    Map.center()
    game.render()

  game.on 'progress', (e) ->
    false

  game.load [
    name: 'island2'
    src: 'resources/map/island2.json'
  ]

window.onresize = ->
  game.resize(window.innerWidth, window.innerHeight)
  Map.center()