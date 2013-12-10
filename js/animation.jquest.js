// Generated by CoffeeScript 1.6.3
var Animation;

Animation = {
  loadChar: function(name, eight) {
    var baseTexture, directions, height, pointX, pointY, textures, width, x, y, _i, _j;
    if (eight == null) {
      eight = true;
    }
    textures = [];
    directions = ['s', 'w', 'e', 'n', 'sw', 'se', 'nw', 'ne'];
    baseTexture = game.cache.getTexture(name).baseTexture;
    width = 32;
    height = 32;
    for (y = _i = 0; _i < 8; y = ++_i) {
      textures[directions[y]] = [];
      for (x = _j = 0; _j < 6; x = ++_j) {
        pointX = x * 32;
        pointY = y * 32;
        textures[directions[y]].push(new gf.Texture(baseTexture, gf.Rectangle(pointX, pointY, width, height)));
      }
    }
    return gf.Sprite({
      n: {
        frames: textures['n']
      },
      s: {
        frames: textures['s']
      },
      w: {
        frames: textures['w']
      },
      e: {
        frames: textures['e']
      },
      nw: {
        frames: textures['nw']
      },
      ne: {
        frames: textures['ne']
      },
      sw: {
        frames: textures['sw']
      },
      se: {
        frames: textures['se']
      }
    }, 1000, 's');
  }
};
