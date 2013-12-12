Input =

  setup: ->
    game.input.keyboard.on gf.Keyboard.KEY.W, Input.up
    game.input.keyboard.on gf.Keyboard.KEY.S, Input.down
    game.input.keyboard.on gf.Keyboard.KEY.A, Input.left
    game.input.keyboard.on gf.Keyboard.KEY.D, Input.right

  up: (e) -> Input.doMove 'n', e.down
  down: (e) -> Input.doMove 's', e.down
  left: (e) -> Input.doMove 'w', e.down
  right: (e) -> Input.doMove 'e', e.down

  doMove: (direction, keyDown) ->
    c = Characters.store['player']
    if keyDown
      switch direction
        when "n" then c.setVelocity 0,c.movespeed
        when "s" then c.setVelocity 0,-c.movespeed
        when "e" then c.setVelocity c.movespeed,0
        when "w" then c.setVelocity -c.movespeed,0
      c.goto(1, direction).play() if c.direction is false
      c.direction = direction

    else
      c.direction = false
      c.goto(1, direction).stop()