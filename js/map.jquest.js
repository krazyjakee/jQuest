var Map;

Map = {
  center: function() {
    var x, y;
    x = Math.floor(world.realSize.x / 2);
    y = Math.floor(world.realSize.y / 2);
    return game.camera.unconstrain().focus(x, y);
  }
};
