// Generated by CoffeeScript 1.6.3
var game, world;

world = game = false;

$(window).ready(function() {
  game = new gf.Game('game', {
    width: window.innerWidth,
    height: window.innerHeight,
    background: 0x000000,
    interactive: true
  });
  game.on('complete', function() {
    game.loadWorld('island2');
    world = game.world;
    Map.center();
    return game.render();
  });
  game.on('progress', function(e) {
    return false;
  });
  return game.load([
    {
      name: 'island2',
      src: 'resources/map/island2.json'
    }
  ]);
});

window.onresize = function() {
  game.resize(window.innerWidth, window.innerHeight);
  return Map.center();
};
