var utils = require('../utils/utils'),
    inherit = require('../utils/inherit'),
    math = require('../math/math'),
    Texture = require('../display/Texture'),
    Vector = require('../math/Vector'),
    PIXI = require('pixi.js');

/**
 * This object represents a tileset used by a Tilemap.
 * There can be multiple Tilesets in a map
 *
 * @class Tileset
 * @extends Texture
 * @constructor
 * @param texture {Texture} The texture to use for the tileset
 * @param settings {Object} All the settings for the tileset
 * @param settings.tilewidth {Number} The width of a single tile in the set
 * @param settings.tileheight {Number} The height of a single tile in the set
 * @param [settings.firstgid=1] {Number} The id of the first tile in the set, defaults to 1
 * @param [settings.spacing=0] {Number} The spacing around tiles in the tileset (in pixels)
 * @param [settings.margin=0] {Number} The margin around a tile in the tileset (in pixels)
 * @param [settings.tileoffset] {Object} The offset to apply to a tile rendered from this tileset
 * @param [settings.tileoffset.x=0] {Number} The X offset to apply to the tile
 * @param [settings.tileoffset.y=0] {Number} The Y offset to apply to the tile
 * @param [settings.properties] {Object} User-defined, custom properties that apply to the tileset
 * @param [settings.tileproperties] {Object} User-defined, custom properties that apply to tiles in the tileset.
 *          The keys of this object should the tile id of the properties
 * @param [settings.imagewidth] {Number} An override for the image width
 * @param [settings.imageheight] {Number} An override for the image height
 */
//TODO: Support external tilesets (TSX files) via the "source" attribute
//see: https://github.com/bjorn/tiled/wiki/TMX-Map-Format#tileset
var Tileset = function(texture, settings) {
    //initialize the base Texture class
    if(texture instanceof Array) {
        this.multi = true;
        Texture.call(this, texture[0].baseTexture);
    } else {
        Texture.call(this, texture.baseTexture || texture);
    }

    //Tiled Editor properties

    /**
     * The first tileId in the tileset
     *
     * @property firstgid
     * @type Number
     */
    this.firstgid = settings.firstgid || 1;

    /**
     * The name of the tileset
     *
     * @property name
     * @type String
     */
    this.name = settings.name;

    /**
     * The size of a tile in the tileset
     *
     * @property tileSize
     * @type Vector
     */
    this.tileSize = new Vector(settings.tilewidth, settings.tileheight);

    /**
     * The spacing around a tile in the tileset
     *
     * @property spacing
     * @type Number
     */
    this.spacing = settings.spacing || 0;

    /**
     * The margin around a tile in the tileset
     *
     * @property margin
     * @type Number
     */
    this.margin = settings.margin || 0;

    /**
     * The offset of tile positions when rendered
     *
     * @property tileoffset
     * @type Number
     */
    this.tileoffset = new Vector(
        settings.tileoffset ? settings.tileoffset.x : 0,
        settings.tileoffset ? settings.tileoffset.y : 0
    );

    //TODO: Support for "terraintypes," "image"
    //see: https://github.com/bjorn/tiled/wiki/TMX-Map-Format#tileset

    //Custom/Optional properties

    /**
     * The number of tiles calculated based on size, margin, and spacing
     *
     * @property numTiles
     * @type Vector
     */
    this.numTiles = this.multi ? texture.length : new Vector(
        math.floor((this.baseTexture.source.width - this.margin) / (this.tileSize.x - this.spacing)),
        math.floor((this.baseTexture.source.height - this.margin) / (this.tileSize.y - this.spacing))
    );

    /**
     * The last tileId in the tileset
     *
     * @property lastgid
     * @type Number
     */
    this.lastgid = this.firstgid + (this.multi ? texture.length : ((this.numTiles.x * this.numTiles.y) || 1)) - 1;

    /**
     * The properties of the tileset
     *
     * @property properties
     * @type Object
     */
    this.properties = settings.properties || {};

    /**
     * The properties of the tiles in the tileset (like collision stuff)
     *
     * @property tileproperties
     * @type Object
     */
    this.tileproperties = settings.tileproperties || {};

    /**
     * The size of the tileset
     *
     * @property size
     * @type Vector
     */
    this.size = this.multi ? Vector.ZERO : new Vector(
        settings.imagewidth || this.baseTexture.source.width,
        settings.imageheight || this.baseTexture.source.height
    );

    /**
     * The texture instances for each tile in the set
     *
     * @property textures
     * @type Array
     */
    this.textures = this.multi ? texture : [];

    //massages strings into the values they should be
    //i.e. "true" becomes the value: true
    this.properties = utils.parseTiledProperties(this.properties);

    //massage tile properties
    for(var k in this.tileproperties) {
        this.tileproperties[k] = utils.parseTiledProperties(this.tileproperties[k]);
    }

    //generate tile textures
    if(!this.multi) {
        for(var t = 0, tl = this.lastgid - this.firstgid + 1; t < tl; ++t) {
            //convert the tileId to x,y coords of the tile in the Texture
            var y = math.floor(t / this.numTiles.x),
                x = (t - (y * this.numTiles.x));

            //get location in pixels
            x = (x * this.tileSize.x) + (x * this.spacing) + this.margin;
            y = (y * this.tileSize.y) + (y * this.spacing) + this.margin;

            this.textures.push(
                new Texture(
                    this.baseTexture,
                    new PIXI.Rectangle(x, y, this.tileSize.x, this.tileSize.y)
                )
            );
        }
    }
};

inherit(Tileset, Texture, {
    /**
     * Gets the tile properties for a tile based on it's ID
     *
     * @method getTileProperties
     * @param tileId {Number} The id of the tile to get the properties for
     * @return {Object} The properties of the tile
     */
    getTileProperties: function(tileId) {
        if(!tileId) return null;

        var flags = Tileset.FLAGS,
            flippedX = tileId & flags.FlippedX,
            flippedY = tileId & flags.FlippedY,
            rotatedCW = tileId & flags.RotatedCW;

        tileId = (tileId & ~Tileset.FLAGS.ALL) - this.firstgid;

        //if less than 0, then this id isn't in this tileset
        if(tileId < 0) return null;

        var props = this.tileproperties[tileId] ?
        //get this value
        this.tileproperties[tileId] :
        //set this id to default values and cache
        this.tileproperties[tileId] = {
            collidable: false,
            breakable: false
        };

        props.flippedX = flippedX;
        props.flippedY = flippedY;
        props.rotatedCW = rotatedCW;

        return props;
    },
    /**
     * Gets the tile texture for a tile based on it's ID
     *
     * @method getTileTexture
     * @param tileId {Number} The id of the tile to get the texture for
     * @return {Texture} The texture for the tile
     */
    getTileTexture: function(tileId) {
        if(!tileId) return null;

        //get the internal ID of the tile in this set (0 indexed)
        tileId = (tileId & ~Tileset.FLAGS.ALL) - this.firstgid;

        //if less than 0, then this id isn't in this tileset
        if(tileId < 0) return null;

        return this.textures[tileId];
    },
    /**
     * Returns whether or not this tileset contains the given tile guid
     *
     * @method contains
     * @param tileId {Number} The ID of the tile to check
     * @return {Boolean}
     */
    contains: function(tileId) {
        if(!tileId) return false;

        tileId &= ~Tileset.FLAGS.ALL;

        return (tileId >= this.firstgid && tileId <= this.lastgid);
    }
});

/**
 * Tileset GID flags, these flags are set on a tile's ID to give it a special property
 *
 * @property FLAGS
 * @static
 */
Tileset.FLAGS = {
    FlippedX: 0x80000000,
    FlippedY: 0x40000000,
    RotatedCW: 0x20000000
};

var mask = 0;
for(var f in Tileset.FLAGS) {
    mask |= Tileset.FLAGS[f];
}

Tileset.FLAGS.ALL = mask;

module.exports = Tileset;
