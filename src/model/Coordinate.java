package model;

public class Coordinate {
	public float x;
	public float y;
	
	public Coordinate(float x, float y) {
		this.x = x;
		this.y = y;
	}
	
	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append("[" + x + ", " + y + "]");
		return sb.toString();
	}
}
