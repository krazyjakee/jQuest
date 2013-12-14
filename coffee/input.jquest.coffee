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

    setvel = (x,y) ->
      v = new gf.Vector(x, y)
      c.velocity = v
      c.setVelocity(v)

    if keyDown
      c.keysDown[direction] = true
      m = c.movespeed
      for k,v of c.keysDown
        if v is true and k isnt direction[0]
          if k is 'n' then direction = 'n' + direction[0]
          if k is 's' then direction = 's' + direction[0]
          if k is 'w' then direction = direction[0] + 'w'
          if k is 'e' then direction = direction[0] + 'e'
      switch direction
        when "n" then setvel(0, -m)
        when "s" then setvel(0, m)
        when "e" then setvel(m, 0)
        when "w" then setvel(-m, 0)
        when 'nw' then setvel(-m, -m)
        when 'ne' then setvel(m, -m)
        when 'sw' then setvel(-m, m)
        when 'se' then setvel(m, m)

      c.goto(1, direction).play() if c.direction isnt direction
      c.direction = direction

    else
      switch direction
        when "n" then setvel(c.velocity.x, 0)
        when "s" then setvel(c.velocity.x, 0)
        when "e" then setvel(0, c.velocity.y)
        when "w" then setvel(0, c.velocity.y)

      c.keysDown[direction] = false

      Keys = false
      for k,v of c.keysDown 
        if v
          Keys = k
          break

      if Keys
        c.goto(1, Keys).play()
      else
        c.goto(1, direction).stop()
      c.direction = false