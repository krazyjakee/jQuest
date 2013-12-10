Input =

  setup: ->
    game.input.keyboard.on gf.Keyboard.KEY.W, Input.up
    game.input.keyboard.on gf.Keyboard.KEY.S, Input.down
    game.input.keyboard.on gf.Keyboard.KEY.A, Input.left
    game.input.keyboard.on gf.Keyboard.KEY.D, Input.right

  up: (e) ->
    false