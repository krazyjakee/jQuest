game = map = tileset = layer = player = cursors = false

preload = ->
  game.load.tilemap 'map', 'resources/map/island2.json', null, Phaser.Tilemap.TILED_JSON
  game.load.image 'tiles', 'resources/img/free_tileset_version_10.png', 32, 32
  game.load.spritesheet 'character', 'resources/img/6Actor_5.png', 32, 32

create = ->

  game.stage.backgroundColor = '#000000'
  map = game.add.tilemap('map')
  map.addTilesetImage('free_tileset_version_10', 'tiles')
  map.setCollisionBetween(1, 13);
  layer = map.createLayer('Ground')
  layer.resizeWorld()

  player = game.add.sprite(4 * 32, 6 * 32, 'character')
  game.camera.follow(player)
  cursors = game.input.keyboard.createCursorKeys()

update = ->
  return

game = new Phaser.Game 800, 600, Phaser.AUTO, 'phaser-game', 
  preload: preload 
  create: create 
  update: update